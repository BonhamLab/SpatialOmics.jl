# SpatialIO/src/spatialdata/from_spatialdata.jl
# Read a SpatialData (Python) OME-ZARR store into a SpatialDataset.
#
# Two strategies:
#   1. Native Zarr.jl (no Python dependency, limited metadata decoding)
#   2. PythonCall bridge (full SpatialData spec, requires spatialdata Python pkg)

using DataFrames
using SparseArrays
using Parquet2

"""
    from_spatialdata(zarr_path::String; use_python::Bool=false) -> SpatialDataset

Read a SpatialData-compatible OME-ZARR store from `zarr_path`.

When `use_python=true`, delegates to the `spatialdata` Python package via
PythonCall.jl for full spec compliance (recommended for complex stores).
When `use_python=false` (default), uses the native Zarr.jl backend.

# Example
```julia
ds = from_spatialdata("/data/my_experiment.zarr")
ds = from_spatialdata("/data/my_experiment.zarr", use_python=true)
```
"""
function from_spatialdata(zarr_path::String; use_python::Bool = false)::SpatialDataset
    if use_python
        return _from_spatialdata_python(zarr_path)
    else
        return _from_spatialdata_native(zarr_path)
    end
end

# ─── Native reader ────────────────────────────────────────────────────────────

function _from_spatialdata_native(zarr_path::String)::SpatialDataset
    isdir(zarr_path) || error("from_spatialdata: path not found: $zarr_path")
    isfile(joinpath(zarr_path, "zarr.json")) ||
        error("from_spatialdata: not a zarr v3 store (no zarr.json): $zarr_path")

    root_meta = _zread_meta(zarr_path)
    ds_meta   = Dict{String,Any}("zarr_attrs" => root_meta)
    ds        = spatial_dataset(; metadata = ds_meta)

    # ── images ───────────────────────────────────────────────────────────────
    images_path = joinpath(zarr_path, "images")
    if isdir(images_path)
        for name in _zlist_groups(images_path)
            img = _read_sd_image(joinpath(images_path, name))
            add_image!(ds, name, img)
        end
    end

    # ── labels ───────────────────────────────────────────────────────────────
    labels_path = joinpath(zarr_path, "labels")
    if isdir(labels_path)
        for name in _zlist_groups(labels_path)
            lbl = _read_sd_labels(joinpath(labels_path, name))
            add_labels!(ds, name, lbl)
        end
    end

    # ── points ───────────────────────────────────────────────────────────────
    points_path = joinpath(zarr_path, "points")
    if isdir(points_path)
        for name in _zlist_groups(points_path)
            pts = _read_sd_points(joinpath(points_path, name))
            add_points!(ds, name, pts)
        end
    end

    # ── shapes ───────────────────────────────────────────────────────────────
    shapes_path = joinpath(zarr_path, "shapes")
    if isdir(shapes_path)
        for name in _zlist_groups(shapes_path)
            shp = _read_sd_shapes(joinpath(shapes_path, name))
            add_shapes!(ds, name, shp)
        end
    end

    # ── tables ───────────────────────────────────────────────────────────────
    tables_path = joinpath(zarr_path, "tables")
    if isdir(tables_path)
        for name in _zlist_groups(tables_path)
            tbl = _read_sd_table(joinpath(tables_path, name))
            add_table!(ds, name, tbl)
        end
    end

    return ds
end

# ─── Image reader ─────────────────────────────────────────────────────────────
#
# Multiscale NGFF images: the group contains numbered sub-groups (0, 1, 2, …)
# for each pyramid level.  We expose level 0 as a lazy ZarrV3Array and stash
# the full metadata for round-trip writing.

function _read_sd_image(path::String)::SpatialImage
    meta = _zread_meta(path)

    # Level-0 array is the full-resolution data
    level0_path = joinpath(path, "0")
    arr = ZarrV3Array(level0_path)

    # Parse axes from NGFF multiscales metadata
    axes_meta = try
        meta.attributes.ome.multiscales[1].axes
    catch
        nothing
    end
    ax_nt = _axes_namedtuple(axes_meta)

    return SpatialImage(arr, ax_nt, Dict{String,Any}("zarr_attrs" => meta))
end

# ─── Labels reader ────────────────────────────────────────────────────────────

function _read_sd_labels(path::String)::SpatialLabels
    meta        = _zread_meta(path)
    level0_path = joinpath(path, "0")
    arr         = ZarrV3Array(level0_path)
    return SpatialLabels(arr, Dict{String,Any}("zarr_attrs" => meta))
end

# ─── Points reader ────────────────────────────────────────────────────────────
#
# Points are stored as a `points.parquet` file inside the group directory.
# Coordinate columns are defined in the group's zarr.json `axes` attribute.
#
# TODO(perf): this is the dominant load-time bottleneck — a Xenium transcript
# file is ~350 MB and is fully materialised into a DataFrame here.  The fix is
# to not force a DataFrame until the user asks for one:
#
#   1. In SpatialOmicsBase, make `features` a parametric type `F` and enforce
#      the Tables.jl Holy-trait in the inner constructor:
#
#        struct SpatialPoints{T<:AbstractFloat, F} <: SpatialElement
#            coordinates::AbstractMatrix{T}
#            features::F
#            metadata::Dict{String,Any}
#            function SpatialPoints(coords, features, meta)
#                Tables.istable(features) ||
#                    throw(ArgumentError("features must satisfy Tables.istable"))
#                new{eltype(coords), typeof(features)}(coords, features, meta)
#            end
#        end
#
#   2. Pass a lazy handle instead of a materialised DataFrame:
#      Arrow.Table (mmap=true) gives zero-copy column access and satisfies
#      Tables.istable.  A thin wrapper around Parquet2.Dataset that reads
#      row-groups on demand would also work.
#
#   3. Widen `coordinates` to `AbstractMatrix{T}` (done in step 1 above) so
#      we can pass a StructArrays view or DiskArray without copying.
#
#   With these changes `from_spatialdata` reads only zarr.json files and
#   returns in <100 ms; data is pulled when the user first indexes features.

function _read_sd_points(path::String)::SpatialPoints
    meta = _zread_meta(path)

    parquet_path = joinpath(path, "points.parquet")
    (isfile(parquet_path) || isdir(parquet_path)) ||
        error("_read_sd_points: missing $parquet_path")

    df = _read_parquet(parquet_path)

    # Determine coordinate columns from metadata (default: x, y, z if present)
    axes = try String.(meta.attributes.axes) catch; ["x", "y"] end
    coord_cols = [a for a in axes if a in names(df)]
    isempty(coord_cols) && (coord_cols = ["x", "y"])

    # Build N×D Float32 coordinate matrix  (hcat of column vectors → N×D)
    coords = Matrix{Float32}(reduce(hcat, [Float32.(df[!, c]) for c in coord_cols]))

    # Feature columns = everything except coordinates
    feat_df = df[:, setdiff(names(df), coord_cols)]

    return SpatialPoints(
        coords,
        feat_df,
        Dict{String,Any}("zarr_attrs" => meta, "coord_cols" => coord_cols),
    )
end

# ─── Shapes reader ────────────────────────────────────────────────────────────
#
# Shapes are stored as a `shapes.parquet` file containing a `geometry` column
# (WKB-encoded polygons/circles) plus optional attribute columns.
#
# TODO(perf): same eager-parquet issue as points — see note above.

function _read_sd_shapes(path::String)::SpatialShapes
    meta = _zread_meta(path)

    parquet_path = joinpath(path, "shapes.parquet")
    (isfile(parquet_path) || isdir(parquet_path)) ||
        error("_read_sd_shapes: missing $parquet_path")

    df = _read_parquet(parquet_path)

    # Keep geometry column as raw bytes in `geometries` and the rest in features
    geom_col = "geometry"
    geoms     = geom_col in names(df) ? df[:, geom_col] : Any[]
    feat_df   = df[:, setdiff(names(df), [geom_col])]

    return SpatialShapes(
        convert(Vector{Any}, geoms),
        feat_df,
        Dict{String,Any}("zarr_attrs" => meta),
    )
end

# ─── Table (AnnData) reader ───────────────────────────────────────────────────
#
# AnnData-style zarr layout:
#   X/           — CSR sparse matrix (data, indices, indptr)
#   obs/         — dataframe (column-order in zarr.json attributes)
#   var/         — dataframe
#   obsm/spatial — 2-D coordinate array
#   uns/         — unstructured metadata

function _read_sd_table(path::String)::SpatialTable
    meta = _zread_meta(path)

    obs = _read_anndata_dataframe(joinpath(path, "obs"))
    var = _read_anndata_dataframe(joinpath(path, "var"))
    X   = _read_csr_matrix(joinpath(path, "X"))

    return SpatialTable(X, obs, var, Dict{String,Any}("zarr_attrs" => meta))
end

# ─── AnnData dataframe reader ─────────────────────────────────────────────────

function _read_anndata_dataframe(path::String)::DataFrame
    isdir(path) || return DataFrame()
    meta = _zread_meta(path)

    col_order = try String.(meta.attributes["column-order"]) catch; String[] end
    index_col = try String(meta.attributes["_index"]) catch; "_index" end

    df = DataFrame()

    # Read index column first (becomes a regular column named "_index")
    idx_path = joinpath(path, index_col)
    if isdir(idx_path) && isfile(joinpath(idx_path, "zarr.json"))
        df[!, "_index"] = _zread_column(idx_path)
    end

    for col in col_order
        col == index_col && continue  # already read
        col_path = joinpath(path, col)
        isdir(col_path) || continue
        isfile(joinpath(col_path, "zarr.json")) || continue
        df[!, col] = _zread_column(col_path)
    end

    return df
end

"""
Read a single column from an AnnData obs/var group.
Handles numeric arrays, string arrays (encoding-type string-array),
and categorical columns (group with categories/ + codes/ sub-arrays).
"""
function _zread_column(col_path::String)
    col_meta = _zread_meta(col_path)

    if String(col_meta.node_type) == "array"
        return _zread_array_1d(col_path)

    elseif String(col_meta.node_type) == "group"
        # Categorical: has categories/ and codes/ sub-arrays
        cats_path  = joinpath(col_path, "categories")
        codes_path = joinpath(col_path, "codes")
        if isdir(cats_path) && isdir(codes_path)
            cats  = _zread_array_1d(cats_path)   # Vector{String}
            codes = _zread_array_1d(codes_path)  # Vector{Int8} (0-based)
            return [cats[Int(c)+1] for c in codes]
        end
        return nothing
    end

    return nothing
end

# ─── CSR sparse matrix reader ─────────────────────────────────────────────────

function _read_csr_matrix(path::String)::SparseMatrixCSC
    isdir(path) || error("_read_csr_matrix: missing group $path")
    meta = _zread_meta(path)

    shape = try Int.(meta.attributes.shape) catch; error("_read_csr_matrix: no shape in $path") end
    n_obs, n_genes = shape[1], shape[2]

    data_vals = Float32.(_zread_array_1d(joinpath(path, "data")))
    indices   = Int.(_zread_array_1d(joinpath(path, "indices")))  # 0-based col idx
    indptr    = Int.(_zread_array_1d(joinpath(path, "indptr")))   # 0-based row ptrs

    # Build COO from CSR then let Julia construct the CSC
    n_rows = n_obs
    I = Vector{Int}(undef, length(data_vals))
    J = Int.(indices) .+ 1   # 0→1-based column indices
    pos = 1
    for row in 1:n_rows
        row_start = indptr[row] + 1
        row_end   = indptr[row+1]
        cnt       = row_end - row_start + 1
        if cnt > 0
            I[pos:pos+cnt-1] .= row
            pos += cnt
        end
    end

    return sparse(I, J, data_vals, n_obs, n_genes)
end

# ─── Helpers ──────────────────────────────────────────────────────────────────

"""
Read a parquet dataset from `path`.  Handles both single `.parquet` files
and partitioned directories (containing `part.N.parquet` files).
"""
function _read_parquet(path::String)::DataFrame
    if isfile(path)
        return DataFrame(Parquet2.Dataset(path))
    elseif isdir(path)
        parts = sort(filter(f -> endswith(f, ".parquet"), readdir(path, join=true)))
        isempty(parts) && return DataFrame()
        return vcat([DataFrame(Parquet2.Dataset(p)) for p in parts]...)
    else
        error("_read_parquet: path not found: $path")
    end
end

function _axes_namedtuple(axes_meta)
    axes_meta === nothing && return NamedTuple()
    names_vec  = [Symbol(String(a.name)) for a in axes_meta]
    units_vec  = [try String(a.unit) catch; "" end for a in axes_meta]
    vals       = Tuple(units_vec)
    return NamedTuple{Tuple(names_vec)}(vals)
end

# ─── Python bridge (stub) ─────────────────────────────────────────────────────

function _from_spatialdata_python(zarr_path::String)::SpatialDataset
    # TODO Phase 1: implement PythonCall bridge
    # sd = pyimport("spatialdata")
    # sdata = sd.read_zarr(zarr_path)
    # ... convert Python objects to Julia SpatialDataset
    error("from_spatialdata (python): not yet implemented")
end
