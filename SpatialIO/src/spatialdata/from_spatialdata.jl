# SpatialIO/src/spatialdata/from_spatialdata.jl
# Read a SpatialData (Python) OME-ZARR store into a SpatialDataset.
#
# Two strategies:
#   1. Native Zarr.jl (no Python dependency, limited metadata decoding)
#   2. PythonCall bridge (full SpatialData spec, requires spatialdata Python pkg)

using DataFrames
using SparseArrays
using Parquet2
using GeometryBasics

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
            ds[name] = img
        end
    end

    # ── labels ───────────────────────────────────────────────────────────────
    labels_path = joinpath(zarr_path, "labels")
    if isdir(labels_path)
        for name in _zlist_groups(labels_path)
            lbl = _read_sd_labels(joinpath(labels_path, name))
            ds[name] = lbl
        end
    end

    # ── points ───────────────────────────────────────────────────────────────
    points_path = joinpath(zarr_path, "points")
    if isdir(points_path)
        for name in _zlist_groups(points_path)
            pts = _read_sd_points(joinpath(points_path, name))
            ds[name] = pts
        end
    end

    # ── shapes ───────────────────────────────────────────────────────────────
    shapes_path = joinpath(zarr_path, "shapes")
    if isdir(shapes_path)
        for name in _zlist_groups(shapes_path)
            shp = _read_sd_shapes(joinpath(shapes_path, name))
            ds[name] = shp
        end
    end

    # ── tables ───────────────────────────────────────────────────────────────
    tables_path = joinpath(zarr_path, "tables")
    if isdir(tables_path)
        for name in _zlist_groups(tables_path)
            tbl = _read_sd_table(joinpath(tables_path, name))
            ds[name] = tbl
        end
    end

    return ds
end

# ─── Image reader ─────────────────────────────────────────────────────────────
#
# Multiscale NGFF images: the group contains numbered sub-groups (0, 1, 2, …)
# for each pyramid level.  We expose level 0 as a lazy ZarrV3Array and stash
# the full metadata for round-trip writing.
#
# If the top-level coordinateTransformations entry is an affine that involves
# an axis permutation or reflection (common for H&E images in Xenium SpatialData
# stores), we detect the normalisation parameters and store them in image
# metadata so that SpatialViz can render the image in the correct global
# coordinate space for overlay with other elements.

function _read_sd_image(path::String)::SpatialImage
    meta = _zread_meta(path)

    # Level-0 array is the full-resolution data
    level0_path = joinpath(path, "0")
    arr = ZarrV3Array(level0_path)

    # Collect additional coarser pyramid levels (1, 2, 3, …)
    pyramid = Any[]
    level_idx = 1
    while isdir(joinpath(path, string(level_idx)))
        push!(pyramid, ZarrV3Array(joinpath(path, string(level_idx))))
        level_idx += 1
    end

    # Parse axes from NGFF multiscales metadata
    axes_meta = if haskey(meta, :attributes) &&
                   haskey(meta.attributes, :ome) &&
                   haskey(meta.attributes.ome, :multiscales) &&
                   !isempty(meta.attributes.ome.multiscales) &&
                   haskey(meta.attributes.ome.multiscales[1], :axes)
        meta.attributes.ome.multiscales[1].axes
    else
        nothing
    end
    ax_nt = _axes_namedtuple(axes_meta)

    img_meta = Dict{String,Any}("zarr_attrs" => meta)

    # Parse coordinate transform; store normalisation parameters in metadata.
    _apply_ngff_transform!(img_meta, meta, size(arr))

    # Fallback: if _apply_ngff_transform! found no transform, check for explicit
    # x_range/y_range attributes written by _write_sd_image_tiles.
    if !haskey(img_meta, "x_range") &&
       haskey(meta, :attributes) &&
       haskey(meta.attributes, :x_range) && haskey(meta.attributes, :y_range)
        xr = meta.attributes.x_range
        yr = meta.attributes.y_range
        img_meta["x_range"] = (Float64(xr[1]), Float64(xr[2]))
        img_meta["y_range"] = (Float64(yr[1]), Float64(yr[2]))
    end

    return SpatialImage(arr, pyramid, ax_nt, img_meta)
end

# ─── NGFF coordinate-transform helpers ────────────────────────────────────────

# Walk the standard NGFF multiscales metadata path (ome.multiscales or multiscales).
function _get_multiscales_meta(meta)
    haskey(meta, :attributes) || return nothing
    a = meta.attributes
    if haskey(a, :ome) && haskey(a.ome, :multiscales) && !isempty(a.ome.multiscales)
        return a.ome.multiscales[1]
    elseif haskey(a, :multiscales) && !isempty(a.multiscales)
        return a.multiscales[1]
    end
    return nothing
end

# Parse the first top-level coordinateTransformations entry.
# Returns (type, data, in_axes, out_axes) or nothing.
function _parse_top_transform(ms)
    haskey(ms, :coordinateTransformations) || return nothing
    ct = ms.coordinateTransformations
    isempty(ct) && return nothing
    t      = ct[1]
    ttype  = String(t.type)
    in_ax  = haskey(t, :input)  ? [String(a.name) for a in t.input.axes]  : nothing
    out_ax = haskey(t, :output) ? [String(a.name) for a in t.output.axes] : nothing
    if ttype == "identity"
        return ("identity", nothing, in_ax, out_ax)
    elseif ttype == "scale"
        return ("scale", Float64.(t.scale), in_ax, out_ax)
    elseif ttype == "affine"
        n_rows = length(t.affine)
        n_cols = length(t.affine[1])
        A = Matrix{Float64}(undef, n_rows, n_cols)
        for i in 1:n_rows, j in 1:n_cols
            A[i, j] = Float64(t.affine[i][j])
        end
        return ("affine", A, in_ax, out_ax)
    else
        return nothing
    end
end

# Inspect the image's NGFF coordinateTransformations and, if a non-trivial
# affine transform is found, write normalisation parameters into `img_meta`:
#   "x_range", "y_range"  — global extent as (Float64, Float64)
#   "perm_yx"             — Bool: local y maps to global x (axes swapped)
#   "flip_dim2"           — Bool: flip post-permutation y dimension
#   "flip_dim3"           — Bool: flip post-permutation x dimension
#
# The +1 offset on global coordinates converts 0-based pixel-centre values
# (used by the NGFF affine) to 1-based Makie pixel coordinates so that
# identity-transform images (rendered at 1..nx) align with affine-transform
# images in the same Makie axis.
function _apply_ngff_transform!(img_meta::Dict{String,Any}, zarr_meta, arr_shape::Tuple)
    ndims_arr = length(arr_shape)
    ndims_arr >= 3 || return   # expect at least (c, ny, nx)

    ms = _get_multiscales_meta(zarr_meta)
    ms === nothing && return

    result = _parse_top_transform(ms)
    result === nothing && return

    ttype, tdata, in_axes, out_axes = result
    ttype == "identity" && return

    if ttype == "scale"
        # Pure scale: physical extents differ from pixel coords.
        # Store extents using the scale factor: pixel i (1-based) is at s*(i-1)+1.
        in_axes === nothing && return
        y_dim = findfirst(==("y"), in_axes)
        x_dim = findfirst(==("x"), in_axes)
        (y_dim === nothing || x_dim === nothing) && return
        ny_l = Float64(arr_shape[end-1])
        nx_l = Float64(arr_shape[end])
        sy = Float64(tdata[y_dim]);  sx = Float64(tdata[x_dim])
        img_meta["x_range"]   = (1.0, sx * (nx_l - 1) + 1.0)
        img_meta["y_range"]   = (1.0, sy * (ny_l - 1) + 1.0)
        img_meta["perm_yx"]   = false
        img_meta["flip_dim2"] = false
        img_meta["flip_dim3"] = false
        return
    end

    ttype == "affine" || return
    (in_axes === nothing || out_axes === nothing) && return

    in_y  = findfirst(==("y"), in_axes);  in_x  = findfirst(==("x"), in_axes)
    out_y = findfirst(==("y"), out_axes); out_x = findfirst(==("x"), out_axes)
    (in_y === nothing || in_x === nothing ||
     out_y === nothing || out_x === nothing) && return

    A = tdata   # (n_out × n_in+1): rows=output dims, cols=input dims + translation

    # 2×2 spatial sub-matrix: output y/x rows × input y/x cols.
    Ayy = A[out_y, in_y];  Ayx = A[out_y, in_x];  ty = A[out_y, end]
    Axy = A[out_x, in_y];  Axx = A[out_x, in_x];  tx = A[out_x, end]

    # Does y_global mainly come from x_local? (axis permutation)
    perm_yx   = abs(Ayx) > abs(Ayy)
    # Sign of dominant coefficient → flip direction
    flip_dim2 = perm_yx ? (Ayx < 0) : (Ayy < 0)   # flip post-permutation y
    flip_dim3 = perm_yx ? (Axy < 0) : (Axx < 0)   # flip post-permutation x

    # Compute global extents from the four image corners (0-based pixel centres).
    ny_l = Float64(arr_shape[end-1])
    nx_l = Float64(arr_shape[end])
    yl   = ny_l - 1.0
    xl   = nx_l - 1.0
    corners_x = [Axy*0  + Axx*0  + tx,
                 Axy*yl + Axx*0  + tx,
                 Axy*0  + Axx*xl + tx,
                 Axy*yl + Axx*xl + tx]
    corners_y = [Ayy*0  + Ayx*0  + ty,
                 Ayy*yl + Ayx*0  + ty,
                 Ayy*0  + Ayx*xl + ty,
                 Ayy*yl + Ayx*xl + ty]

    # +1 shifts from 0-based pixel centres to 1-based Makie pixel coordinates.
    img_meta["x_range"]   = (minimum(corners_x) + 1.0, maximum(corners_x) + 1.0)
    img_meta["y_range"]   = (minimum(corners_y) + 1.0, maximum(corners_y) + 1.0)
    img_meta["perm_yx"]   = perm_yx
    img_meta["flip_dim2"] = flip_dim2
    img_meta["flip_dim3"] = flip_dim3
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
# Coordinates are read eagerly (small); all other feature columns are kept as
# a lazy Parquet2.Dataset that is materialised only when the user accesses a
# spatial view, crop, or writes back to disk.

function _read_sd_points(path::String)::SpatialPoints
    meta = _zread_meta(path)

    parquet_path = joinpath(path, "points.parquet")
    (isfile(parquet_path) || isdir(parquet_path)) ||
        error("_read_sd_points: missing $parquet_path")

    # Read all rows eagerly. Parquet2.Dataset(dir) returns 0 rows via Tables.jl
    # for partitioned datasets, so we use _read_parquet (part-by-part vcat).
    features_df = _read_parquet(parquet_path)

    # Determine coordinate columns from metadata (default: x, y)
    all_cols = names(features_df)
    axes = haskey(meta, :attributes) && haskey(meta.attributes, :axes) ?
           String.(meta.attributes.axes) : ["x", "y"]
    coord_cols = [a for a in axes if a in all_cols]
    isempty(coord_cols) && (coord_cols = ["x", "y"])

    # Build N×D Float32 coordinate matrix
    coords = Matrix{Float32}(reduce(hcat, [Float32.(features_df[!, c]) for c in coord_cols]))

    # Apply the element-level coordinate transform (local → global).
    n_coord = length(coord_cols)
    scales  = _parse_element_scale(meta, n_coord)
    if any(!=(1.0), scales)
        for (i, s) in enumerate(scales)
            s == 1.0 && continue
            coords[:, i] .*= Float32(s)
        end
    end

    return SpatialPoints(
        coords,
        features_df,
        Dict{String,Any}("zarr_attrs" => meta, "coord_cols" => coord_cols),
    )
end

# ─── Element-level coordinate transform helper ───────────────────────────────
#
# SpatialData stores a scale (or identity) coordinateTransformation in each
# element's zarr.json attributes.  These map from the element's local coordinate
# system to a named global coordinate system so that all elements align when
# displayed together.
#
# Returns (sx, sy) scale factors.  Identity / missing transform → (1.0, 1.0).
# For points the third axis (z) is also returned as sz.

function _parse_element_scale(meta, n_spatial::Int=2)::Vector{Float64}
    ones_out = ones(Float64, n_spatial)
    haskey(meta, :attributes) || return ones_out
    attrs = meta.attributes
    haskey(attrs, :coordinateTransformations) || return ones_out
    ct = attrs.coordinateTransformations
    isempty(ct) && return ones_out
    t = ct[1]
    String(t.type) == "scale" || return ones_out
    raw = Float64.(t.scale)
    # Return the first n_spatial values.
    length(raw) >= n_spatial || return ones_out
    return raw[1:n_spatial]
end

# Scale a single GeometryBasics geometry by (sx, sy).
_apply_scale(::Nothing, ::Float64, ::Float64) = nothing
_apply_scale(g, ::Float64, ::Float64)         = g   # unknown type — pass through

function _apply_scale(g::GeometryBasics.Polygon, sx::Float64, sy::Float64)
    scale_pt(p) = GeometryBasics.Point2f(Float32(p[1] * sx), Float32(p[2] * sy))
    ext  = scale_pt.(g.exterior)
    ints = [scale_pt.(ring) for ring in g.interiors]
    return GeometryBasics.Polygon(ext, ints)
end

function _apply_scale(g::GeometryBasics.Circle, sx::Float64, sy::Float64)
    r_scale = (sx + sy) / 2.0   # uniform when sx ≈ sy; graceful otherwise
    return GeometryBasics.Circle(
        GeometryBasics.Point2f(Float32(g.center[1] * sx), Float32(g.center[2] * sy)),
        Float32(g.r * r_scale),
    )
end

# ─── Shapes reader ────────────────────────────────────────────────────────────
#
# Shapes are stored as a `shapes.parquet` file with a `geometry` WKB column
# plus optional attribute columns.
#
# SpatialData uses two geometry conventions:
#   • Polygon (WKB type 3) — cell boundaries, tissue regions
#   • Point + radius column  — circles (cell spots, Visium spots, "cell_circles")
#     The WKB encodes only the centre; the circle is reconstructed as a polygon.
#
# WKB geometry bytes are decoded inline to GeometryBasics.Polygon so that
# SpatialViz can render them directly with poly!.

function _read_sd_shapes(path::String)::SpatialShapes
    meta = _zread_meta(path)

    parquet_path = joinpath(path, "shapes.parquet")
    (isfile(parquet_path) || isdir(parquet_path)) ||
        error("_read_sd_shapes: missing $parquet_path")

    df = _read_parquet(parquet_path)

    geom_col   = "geometry"
    has_radius = "radius" in names(df)

    geoms = if geom_col in names(df)
        if has_radius
            # Circle format: WKB Point centre + separate radius column
            Any[_decode_wkb_circle(bytes, r)
                for (bytes, r) in zip(df[:, geom_col], df.radius)]
        else
            Any[_decode_wkb(bytes) for bytes in df[:, geom_col]]
        end
    else
        Any[]
    end

    # Apply the element-level coordinate transform (local → global) so that
    # shapes align with images that carry the same global coordinate system.
    sx, sy = _parse_element_scale(meta, 2)
    if sx != 1.0 || sy != 1.0
        geoms = Any[_apply_scale(g, sx, sy) for g in geoms]
    end

    feat_df = df[:, setdiff(names(df), [geom_col])]
    return SpatialShapes(geoms, feat_df, Dict{String,Any}("zarr_attrs" => meta))
end

# ─── WKB decoder ──────────────────────────────────────────────────────────────
#
# Supports Polygon (type 3). Only little-endian byte order (0x01) is handled;
# Z/M flags (0x80000000, 0x40000000) are stripped before type comparison.

function _decode_wkb(bytes::Vector{UInt8})::Union{GeometryBasics.Polygon, Nothing}
    isempty(bytes) && return nothing
    io = IOBuffer(bytes)

    byte_order = Base.read(io, UInt8)
    byte_order == 0x01 || return nothing   # only little-endian supported

    geom_type = ltoh(Base.read(io, UInt32)) & 0x0000FFFF   # strip Z/M/SRID flags
    geom_type == UInt32(3) || return nothing                # we only decode Polygon

    num_rings = Int(ltoh(Base.read(io, UInt32)))
    num_rings == 0 && return nothing

    exterior = _read_wkb_ring(io)
    holes = Vector{Vector{GeometryBasics.Point2f}}()
    for _ in 2:num_rings
        push!(holes, _read_wkb_ring(io))
    end

    isempty(exterior) && return nothing
    return GeometryBasics.Polygon(exterior, holes)
end

# Decode a WKB Point (type 1) centre + radius into a GeometryBasics.Circle.
function _decode_wkb_circle(
    bytes  ::Vector{UInt8},
    radius ::Real,
)::Union{GeometryBasics.Circle, Nothing}
    isempty(bytes) && return nothing
    io = IOBuffer(bytes)

    byte_order = Base.read(io, UInt8)
    byte_order == 0x01 || return nothing

    geom_type = ltoh(Base.read(io, UInt32)) & 0x0000FFFF
    geom_type == UInt32(1) || return nothing   # must be Point

    cx = Float32(ltoh(Base.read(io, Float64)))
    cy = Float32(ltoh(Base.read(io, Float64)))
    return GeometryBasics.Circle(GeometryBasics.Point2f(cx, cy), Float32(radius))
end

function _read_wkb_ring(io::IOBuffer)
    n = Int(ltoh(Base.read(io, UInt32)))
    pts = Vector{GeometryBasics.Point2f}(undef, n)
    for i in 1:n
        x = Float32(ltoh(Base.read(io, Float64)))
        y = Float32(ltoh(Base.read(io, Float64)))
        pts[i] = GeometryBasics.Point2f(x, y)
    end
    return pts
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

    attrs     = haskey(meta, :attributes) ? meta.attributes : nothing
    col_order = attrs !== nothing && haskey(attrs, "column-order") ?
                String.(attrs["column-order"]) : String[]
    index_col = attrs !== nothing && haskey(attrs, "_index") ?
                String(attrs["_index"]) : "_index"

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

    haskey(meta, :attributes) && haskey(meta.attributes, :shape) ||
        error("_read_csr_matrix: no shape in $path")
    shape = Int.(meta.attributes.shape)
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
Open a parquet file/directory as a lazy `Parquet2.Dataset` (satisfies `Tables.istable`).
For partitioned directories, returns the first part's dataset for column-name discovery
and lazy reads; full materialisation via `_read_parquet` is used when all rows are needed.
"""
function _open_parquet(path::String)::Parquet2.Dataset
    if isfile(path)
        return Parquet2.Dataset(path)
    elseif isdir(path)
        parts = sort(filter(f -> endswith(f, ".parquet"), readdir(path, join=true)))
        isempty(parts) && error("_open_parquet: no .parquet files in $path")
        # Return the first part; partitioned read happens in _read_parquet.
        return Parquet2.Dataset(parts[1])
    else
        error("_open_parquet: path not found: $path")
    end
end

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

# Like _read_parquet but reads only the specified columns (avoids loading all columns
# for large datasets where only coordinates are needed).
function _read_parquet_cols(path::String, cols::Vector{Symbol})::DataFrame
    if isfile(path)
        return DataFrame(Parquet2.Dataset(path; readercolumns = cols))
    elseif isdir(path)
        parts = sort(filter(f -> endswith(f, ".parquet"), readdir(path, join=true)))
        isempty(parts) && return DataFrame()
        return vcat([DataFrame(Parquet2.Dataset(p; readercolumns = cols)) for p in parts]...)
    else
        error("_read_parquet_cols: path not found: $path")
    end
end

function _axes_namedtuple(axes_meta)
    axes_meta === nothing && return NamedTuple()
    names_vec  = [Symbol(String(a.name)) for a in axes_meta]
    units_vec  = [haskey(a, :unit) ? String(a.unit) : "" for a in axes_meta]
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
