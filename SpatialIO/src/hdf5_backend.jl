# SpatialIO/src/hdf5_backend.jl
#
# Native HDF5 persistence for SpatialDataset.
#
# On-disk layout
# ──────────────
#   /                            HDF5 root
#     metadata (string attr)     JSON-serialised ds.metadata
#     images/<name>/
#       data                     N-D float array (collected into memory)
#       axes     (string attr)   JSON {"names":[...], "units":[...]}
#       metadata (string attr)   JSON
#     labels/<name>/
#       data                     N-D integer array
#       metadata (string attr)   JSON
#     points/<name>/
#       coordinates              N×D Float32 matrix  (rows = points)
#       features/<col>           1-D array per feature column
#       coord_cols  (string attr) JSON array of column names
#       feat_cols   (string attr) JSON array of feature column names
#       metadata    (string attr) JSON
#     shapes/<name>/
#       features/<col>           1-D array per feature column
#       feat_cols   (string attr) JSON
#       metadata    (string attr) JSON
#       (WKB geometries are not stored in this backend; round-trip via Zarr instead)
#     tables/<name>/
#       X/
#         data                   Float32 NNZ values  (CSC format)
#         indices                Int32 row-indices (0-based)
#         indptr                 Int32 column pointers (0-based)
#         shape   (string attr)  JSON [n_obs, n_genes]
#         format  (string attr)  "csc"
#       obs/<col>                1-D array per obs column
#       obs_cols (string attr)   JSON column order
#       var/<col>                1-D array per var column
#       var_cols (string attr)   JSON column order
#       metadata (string attr)   JSON

using HDF5
using JSON
using SparseArrays

# ─── Public API ───────────────────────────────────────────────────────────────

"""
    write_hdf5(ds::SpatialDataset, path::String)

Serialise `ds` to a native HDF5 file at `path`.

All element arrays are collected into memory before writing (no lazy semantics
on the write path). Use `write_zarr` / `to_spatialdata` for large datasets that
benefit from chunked on-disk storage.

# Example
```julia
write_hdf5(ds, "my_dataset.h5")
```
"""
function write_hdf5(ds::SpatialDataset, path::String)
    HDF5.h5open(path, "w") do fid
        HDF5.write_attribute(fid, "metadata", JSON.json(ds.metadata))

        # ── images ────────────────────────────────────────────────────────────
        if !isempty(ds.images)
            ig = HDF5.create_group(fid, "images")
            for (name, img) in ds.images
                g = HDF5.create_group(ig, name)
                g["data"] = collect(img.data)
                axes_names = [String(k) for k in keys(img.axes)]
                axes_units = [String(v) for v in values(img.axes)]
                HDF5.write_attribute(g, "axes",
                    JSON.json(Dict("names" => axes_names, "units" => axes_units)))
                HDF5.write_attribute(g, "metadata", JSON.json(img.metadata))
            end
        end

        # ── labels ────────────────────────────────────────────────────────────
        if !isempty(ds.labels)
            lg = HDF5.create_group(fid, "labels")
            for (name, lbl) in ds.labels
                g = HDF5.create_group(lg, name)
                g["data"] = collect(lbl.data)
                HDF5.write_attribute(g, "metadata", JSON.json(lbl.metadata))
            end
        end

        # ── points ────────────────────────────────────────────────────────────
        if !isempty(ds.points)
            pg = HDF5.create_group(fid, "points")
            for (name, pts) in ds.points
                g = HDF5.create_group(pg, name)
                g["coordinates"] = Matrix{Float32}(pts.coordinates)
                g["labels"]      = pts.labels   # Vector{String}; empty if no labels
                if pts.features !== nothing
                    fg = HDF5.create_group(g, "features")
                    for (col, vec) in pairs(pts.features)
                        _hdf5_write_column(fg, String(col), vec)
                    end
                end
                HDF5.write_attribute(g, "metadata", JSON.json(pts.metadata))
            end
        end

        # ── shapes ────────────────────────────────────────────────────────────
        if !isempty(ds.shapes)
            sg = HDF5.create_group(fid, "shapes")
            for (name, shp) in ds.shapes
                g      = HDF5.create_group(sg, name)
                fg     = HDF5.create_group(g, "features")
                cols   = Tables.columns(shp)
                fnames = String.(keys(cols))
                for (fname, fvals) in pairs(cols)
                    _hdf5_write_column(fg, String(fname), collect(fvals))
                end
                HDF5.write_attribute(g, "feat_cols", JSON.json(fnames))
                HDF5.write_attribute(g, "metadata",  JSON.json(shp.metadata))
            end
        end

        # ── tables ────────────────────────────────────────────────────────────
        if !isempty(ds.tables)
            tg = HDF5.create_group(fid, "tables")
            for (name, tbl) in ds.tables
                g  = HDF5.create_group(tg, name)
                _hdf5_write_csc(g, tbl.data)
                _hdf5_write_dataframe(g, "obs", tbl.obs, "obs_cols")
                _hdf5_write_dataframe(g, "var", tbl.var, "var_cols")
                HDF5.write_attribute(g, "metadata", JSON.json(tbl.metadata))
            end
        end
    end
end

"""
    read_hdf5(path::String) -> SpatialDataset

Read a `SpatialDataset` from a native HDF5 file written by `write_hdf5`.

# Example
```julia
ds = read_hdf5("my_dataset.h5")
```
"""
function read_hdf5(path::String)::SpatialDataset
    isfile(path) || error("read_hdf5: file not found: $path")

    HDF5.h5open(path, "r") do fid
        meta_str = HDF5.haskey(HDF5.attributes(fid), "metadata") ?
            HDF5.read_attribute(fid, "metadata") : "{}"
        meta = Dict{String,Any}(JSON.parse(meta_str))
        ds   = spatial_dataset(; metadata = meta)

        # ── images ────────────────────────────────────────────────────────────
        if HDF5.haskey(fid, "images")
            for name in keys(fid["images"])
                g    = fid["images"][name]
                data = Base.read(g["data"])
                axes_str = HDF5.haskey(HDF5.attributes(g), "axes") ?
                    HDF5.read_attribute(g, "axes") : "{}"
                axes_d = JSON.parse(axes_str)
                ax_names = get(axes_d, "names", String[])
                ax_units = get(axes_d, "units", String[])
                ax_nt    = _build_axes_nt(ax_names, ax_units)
                img_meta = _read_json_attr(g, "metadata")
                ds[name] = SpatialImage(data, ax_nt, img_meta)
            end
        end

        # ── labels ────────────────────────────────────────────────────────────
        if HDF5.haskey(fid, "labels")
            for name in keys(fid["labels"])
                g    = fid["labels"][name]
                data = Base.read(g["data"])
                lbl_meta = _read_json_attr(g, "metadata")
                ds[name] = SpatialLabels(data, lbl_meta)
            end
        end

        # ── points ────────────────────────────────────────────────────────────
        if HDF5.haskey(fid, "points")
            for name in keys(fid["points"])
                g      = fid["points"][name]
                coords = Matrix{Float32}(Base.read(g["coordinates"]))
                labels = HDF5.haskey(g, "labels") ?
                         Vector{String}(Base.read(g["labels"])) : String[]
                features = nothing
                if HDF5.haskey(g, "features")
                    fg    = g["features"]
                    names = Tuple(Symbol.(keys(fg)))
                    vecs  = Tuple(_hdf5_read_column(fg[String(n)]) for n in names)
                    features = NamedTuple{names}(vecs)
                end
                pts_meta = _read_json_attr(g, "metadata")
                ds[name] = SpatialPoints{Float32}(coords, labels, features, pts_meta)
            end
        end

        # ── shapes ────────────────────────────────────────────────────────────
        if HDF5.haskey(fid, "shapes")
            for name in keys(fid["shapes"])
                g = fid["shapes"][name]
                feat_cols = HDF5.haskey(HDF5.attributes(g), "feat_cols") ?
                    Vector{String}(JSON.parse(HDF5.read_attribute(g, "feat_cols"))) :
                    String[]
                feat_df = DataFrame()
                if HDF5.haskey(g, "features")
                    fg = g["features"]
                    for col in feat_cols
                        HDF5.haskey(fg, col) || continue
                        feat_df[!, col] = _hdf5_read_column(fg[col])
                    end
                end
                shp_meta = _read_json_attr(g, "metadata")
                n = nrow(feat_df)
                ds[name] = SpatialShapes(fill(nothing, n), feat_df, shp_meta)
            end
        end

        # ── tables ────────────────────────────────────────────────────────────
        if HDF5.haskey(fid, "tables")
            for name in keys(fid["tables"])
                g   = fid["tables"][name]
                X   = _hdf5_read_csc(g)
                obs = _hdf5_read_dataframe(g, "obs", "obs_cols")
                var = _hdf5_read_dataframe(g, "var", "var_cols")
                tbl_meta = _read_json_attr(g, "metadata")
                ds[name] = SpatialTable(X, obs, var, tbl_meta)
            end
        end

        return ds
    end
end

# ─── Internal helpers ─────────────────────────────────────────────────────────

function _hdf5_write_column(grp, col::String, vals::AbstractVector)
    if eltype(vals) <: AbstractString
        grp[col] = Vector{String}(vals)
    elseif eltype(vals) <: Real
        grp[col] = collect(vals)
    else
        # Fallback: JSON-encode heterogeneous columns
        grp[col] = JSON.json(collect(vals))
    end
end

function _hdf5_read_column(dset::HDF5.Dataset)
    val = Base.read(dset)
    # JSON-encoded fallback was written as a scalar string
    val isa String && return JSON.parse(val)
    return val
end

function _hdf5_write_csc(grp::HDF5.Group, X::AbstractMatrix)
    xg = HDF5.create_group(grp, "X")
    M  = isa(X, SparseMatrixCSC) ? X : sparse(X)
    xg["data"]    = Vector{Float32}(M.nzval)
    xg["indices"] = Vector{Int32}(M.rowval .- 1)
    xg["indptr"]  = Vector{Int32}(M.colptr .- 1)
    HDF5.write_attribute(xg, "shape",  JSON.json([size(M, 1), size(M, 2)]))
    HDF5.write_attribute(xg, "format", "csc")
end

function _hdf5_read_csc(grp::HDF5.Group)::SparseMatrixCSC
    HDF5.haskey(grp, "X") || return spzeros(Float32, 0, 0)
    xg      = grp["X"]
    vals    = Vector{Float32}(Base.read(xg["data"]))
    rowval  = Vector{Int}(Base.read(xg["indices"])) .+ 1   # 0→1-based
    colptr  = Vector{Int}(Base.read(xg["indptr"]))  .+ 1
    shape   = Vector{Int}(JSON.parse(HDF5.read_attribute(xg, "shape")))
    return SparseMatrixCSC(shape[1], shape[2], colptr, rowval, vals)
end

function _hdf5_write_dataframe(
    grp      ::HDF5.Group,
    subname  ::String,
    df       ::DataFrame,
    cols_attr::String,
)
    g    = HDF5.create_group(grp, subname)
    cols = names(df)
    for col in cols
        _hdf5_write_column(g, col, df[!, col])
    end
    HDF5.write_attribute(grp, cols_attr, JSON.json(cols))
end

function _hdf5_read_dataframe(
    grp      ::HDF5.Group,
    subname  ::String,
    cols_attr::String,
)::DataFrame
    HDF5.haskey(grp, subname) || return DataFrame()
    g    = grp[subname]
    cols = HDF5.haskey(HDF5.attributes(grp), cols_attr) ?
        Vector{String}(JSON.parse(HDF5.read_attribute(grp, cols_attr))) :
        String.(keys(g))
    df   = DataFrame()
    for col in cols
        HDF5.haskey(g, col) || continue
        df[!, col] = _hdf5_read_column(g[col])
    end
    return df
end

function _read_json_attr(obj, key::String)::Dict{String,Any}
    HDF5.haskey(HDF5.attributes(obj), key) || return Dict{String,Any}()
    Dict{String,Any}(JSON.parse(HDF5.read_attribute(obj, key)))
end

function _build_axes_nt(names::Vector, units::Vector)
    isempty(names) && return NamedTuple()
    syms = Tuple(Symbol.(names))
    vals = Tuple(String.(units))
    NamedTuple{syms}(vals)
end
