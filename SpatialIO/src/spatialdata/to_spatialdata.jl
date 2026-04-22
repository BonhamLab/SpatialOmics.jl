# SpatialIO/src/spatialdata/to_spatialdata.jl
# Write a SpatialDataset to a SpatialData-compatible OME-ZARR store.
using DataFrames
using SparseArrays
using Parquet2

"""
    to_spatialdata(ds::SpatialDataset, zarr_path::String; use_python::Bool=false)

Serialise `ds` as a SpatialData-compatible OME-ZARR store at `zarr_path`.

When `use_python=true`, uses the `spatialdata` Python package for writing
(guarantees spec compliance); otherwise uses the native Zarr.jl backend.

# Example
```julia
to_spatialdata(ds, "/output/my_experiment.zarr")
```
"""
function to_spatialdata(
    ds         ::SpatialDataset,
    zarr_path  ::String;
    use_python ::Bool = false,
)
    if use_python
        _to_spatialdata_python(ds, zarr_path)
    else
        _to_spatialdata_native(ds, zarr_path)
    end
    return zarr_path
end

# ─── Helpers ─────────────────────────────────────────────────────────────────

# Extract the group-level attributes dict from an element's stored zarr metadata,
# falling back to an empty dict if absent (elements constructed outside of a Zarr
# round-trip won't have zarr_attrs).
_zarr_group_attrs(meta::Dict{String,Any}) =
    let za = get(meta, "zarr_attrs", nothing)
        za !== nothing && haskey(za, "attributes") ? za["attributes"] : Dict{String,Any}()
    end

# ─── Native writer ────────────────────────────────────────────────────────────

function _to_spatialdata_native(ds::SpatialDataset, zarr_path::String)
    mkpath(zarr_path)

    # Write root zarr.json
    root_attrs = _zarr_group_attrs(ds.metadata)
    _zwrite_group(zarr_path, root_attrs)

    # ── images ────────────────────────────────────────────────────────────────
    if !isempty(ds.images)
        _zwrite_group(joinpath(zarr_path, "images"))
        for (name, img) in ds.images
            _write_sd_image(joinpath(zarr_path, "images", name), img)
        end
    end

    # ── labels ────────────────────────────────────────────────────────────────
    if !isempty(ds.labels)
        _zwrite_group(joinpath(zarr_path, "labels"))
        for (name, lbl) in ds.labels
            _write_sd_labels(joinpath(zarr_path, "labels", name), lbl)
        end
    end

    # ── points ────────────────────────────────────────────────────────────────
    if !isempty(ds.points)
        _zwrite_group(joinpath(zarr_path, "points"))
        for (name, pts) in ds.points
            _write_sd_points(joinpath(zarr_path, "points", name), pts)
        end
    end

    # ── shapes ────────────────────────────────────────────────────────────────
    if !isempty(ds.shapes)
        _zwrite_group(joinpath(zarr_path, "shapes"))
        for (name, shp) in ds.shapes
            _write_sd_shapes(joinpath(zarr_path, "shapes", name), shp)
        end
    end

    # ── tables ────────────────────────────────────────────────────────────────
    if !isempty(ds.tables)
        _zwrite_group(joinpath(zarr_path, "tables"))
        for (name, tbl) in ds.tables
            _write_sd_table(joinpath(zarr_path, "tables", name), tbl)
        end
    end
end

# ─── Image writer ─────────────────────────────────────────────────────────────

function _write_sd_image(path::String, img::SpatialImage)
    mkpath(path)
    group_attrs = _zarr_group_attrs(img.metadata)
    _zwrite_group(path, group_attrs)

    data = img.data
    N    = ndims(data)
    cs   = if data isa ZarrV3Array
        data.chunk_shape
    else
        NTuple{N,Int}(ntuple(i -> i == 1 ? size(data, 1) : min(size(data, i), 4096), N))
    end
    _zwrite_array(joinpath(path, "0"), collect(data); chunk_shape = cs)
end

# ─── Labels writer ────────────────────────────────────────────────────────────

function _write_sd_labels(path::String, lbl::SpatialLabels)
    mkpath(path)
    group_attrs = _zarr_group_attrs(lbl.metadata)
    _zwrite_group(path, group_attrs)

    data = lbl.data
    cs   = data isa ZarrV3Array ? data.chunk_shape : size(data)
    _zwrite_array(joinpath(path, "0"), collect(data); chunk_shape = NTuple{ndims(data),Int}(cs))
end

# ─── Points writer ────────────────────────────────────────────────────────────

function _write_sd_points(path::String, pts::SpatialPoints)
    mkpath(path)
    group_attrs = _zarr_group_attrs(pts.metadata)
    _zwrite_group(path, group_attrs)

    coord_cols = haskey(pts.metadata, "coord_cols") ?
                 String.(pts.metadata["coord_cols"]) : ["x", "y"]
    D = size(pts.coordinates, 2)
    coord_cols = coord_cols[1:min(length(coord_cols), D)]

    # Build DataFrame: coordinates + features (features may be a lazy table)
    df = DataFrame(pts.coordinates, coord_cols)
    for col in Tables.columnnames(pts.features)
        df[!, String(col)] = Tables.getcolumn(pts.features, col)
    end

    Parquet2.writefile(joinpath(path, "points.parquet"), df)
end

# ─── WKB polygon encoder ─────────────────────────────────────────────────────

"""
    _encode_wkb_polygon(geom::GeometryBasics.Polygon) -> Vector{UInt8}

Encode a polygon as little-endian WKB. Format:
  1 byte  byte order mark (0x01 = little-endian)
  4 bytes geometry type   (0x03000000 = Polygon)
  4 bytes ring count      (1 for exterior-only polygons)
  For each ring:
    4 bytes   point count
    n×16 bytes  Float64 x, y pairs

Matches the format decoded by `_decode_wkb` in `from_spatialdata.jl`.
"""
function _encode_wkb_polygon(geom::GeometryBasics.Polygon)
    io = IOBuffer()
    Base.write(io, UInt8(0x01))         # little-endian
    Base.write(io, UInt32(3))           # WKB type: polygon
    ring = GeometryBasics.coordinates(geom)
    Base.write(io, UInt32(1))           # one ring (exterior only)
    Base.write(io, UInt32(length(ring)))
    for pt in ring
        Base.write(io, Float64(pt[1]))
        Base.write(io, Float64(pt[2]))
    end
    return take!(io)
end

# ─── Shapes writer ────────────────────────────────────────────────────────────

function _write_sd_shapes(path::String, shp::SpatialShapes{G,D}) where {G,D}
    mkpath(path)
    group_attrs = _zarr_group_attrs(shp.metadata)
    _zwrite_group(path, group_attrs)

    df = isempty(shp.shapes) ? DataFrame() : DataFrame(shp)

    if !isempty(shp.shapes)
        first_geom = shp.shapes[1].geometry
        if first_geom isa GeometryBasics.Polygon
            df[!, "geometry"] = [_encode_wkb_polygon(s.geometry) for s in shp.shapes]
        elseif first_geom isa GeometryBasics.Circle
            @warn "_write_sd_shapes: Circle WKB encoding not yet implemented — skipping geometry column"
        elseif first_geom !== nothing
            @warn "_write_sd_shapes: unsupported geometry type $(typeof(first_geom)) — skipping geometry column"
        end
    end

    Parquet2.writefile(joinpath(path, "shapes.parquet"), df)
end

# ─── Table (AnnData) writer ───────────────────────────────────────────────────

function _write_sd_table(path::String, tbl::SpatialTable)
    mkpath(path)
    group_attrs = _zarr_group_attrs(tbl.metadata)
    _zwrite_group(path, group_attrs)

    _write_anndata_dataframe(joinpath(path, "obs"), tbl.obs)
    _write_anndata_dataframe(joinpath(path, "var"), tbl.var)
    _write_csr_matrix(joinpath(path, "X"), tbl.data)

    # Write empty required sub-groups
    for sub in ("obsm", "varm", "obsp", "varp", "layers", "uns")
        _zwrite_group(joinpath(path, sub))
    end
    # raw group
    _zwrite_group(joinpath(path, "raw"))
end

# ─── AnnData dataframe writer ─────────────────────────────────────────────────

function _write_anndata_dataframe(path::String, df::DataFrame)
    mkpath(path)
    cols = names(df)

    # Determine which columns are string vs numeric
    col_order = [c for c in cols if c != "_index"]

    obs_attrs = Dict{String,Any}(
        "column-order"    => col_order,
        "_index"          => "_index",
        "encoding-type"   => "dataframe",
        "encoding-version" => "0.2.0",
    )
    _zwrite_group(path, obs_attrs)

    for col in cols
        col_path = joinpath(path, col)
        vals = df[!, col]

        # Strip Union{Missing, T} → T by replacing missing with a sentinel value.
        if Missing <: eltype(vals)
            T_base = nonmissingtype(eltype(vals))
            if T_base <: AbstractFloat
                vals = coalesce.(vals, T_base(NaN))
            elseif T_base <: Integer
                vals = coalesce.(vals, zero(T_base))
            elseif T_base <: AbstractString || T_base == String
                vals = coalesce.(vals, "")
            else
                vals = coalesce.(vals, zero(T_base))
            end
        end

        if eltype(vals) <: AbstractString || eltype(vals) == String
            _zwrite_string_array(col_path, Vector{String}(vals))
            _write_string_array_meta(col_path)
        else
            _zwrite_array(col_path, collect(vals);
                          chunk_shape = (length(vals),))
        end
    end
end

function _write_string_array_meta(path::String)
    # Patch attributes onto existing zarr.json for string-array encoding
    meta_path = joinpath(path, "zarr.json")
    meta = JSON.parse(Base.read(meta_path, String))
    meta["attributes"] = Dict{String,Any}(
        "encoding-type"    => "string-array",
        "encoding-version" => "0.2.0",
    )
    Base.write(meta_path, JSON.json(meta))
end

# ─── CSR matrix writer ────────────────────────────────────────────────────────

function _write_csr_matrix(path::String, X::AbstractMatrix)
    mkpath(path)
    n_obs, n_genes = size(X)

    # Convert SparseMatrixCSC → CSR for writing
    # Julia stores CSC; we need CSR for spatialdata compatibility
    # CSC: A[i,j] col-major → CSR: A[i,j] row-major
    # Transpose trick: CSR(A) = CSC(A')
    Xt = sparse(X')  # SparseMatrixCSC{T, Ti} where Xt is n_genes × n_obs (CSC of A^T = CSR of A)

    data_vals = Float32.(Xt.nzval)
    indices   = Int32.(Xt.rowval .- 1)  # 1-based → 0-based row indices (= col indices of A)
    indptr    = Int32.(Xt.colptr .- 1)  # 1-based → 0-based col ptrs (= row ptrs of CSR of A)

    x_attrs = Dict{String,Any}(
        "shape"            => [n_obs, n_genes],
        "encoding-type"    => "csr_matrix",
        "encoding-version" => "0.1.0",
    )
    _zwrite_group(path, x_attrs)

    _zwrite_array(joinpath(path, "data"),    data_vals; chunk_shape = (length(data_vals),))
    _zwrite_array(joinpath(path, "indices"), indices;   chunk_shape = (length(indices),))
    _zwrite_array(joinpath(path, "indptr"),  indptr;    chunk_shape = (length(indptr),))
end

# ─── Python bridge (stub) ─────────────────────────────────────────────────────

function _to_spatialdata_python(ds::SpatialDataset, zarr_path::String)
    error("to_spatialdata (python): not yet implemented")
end
