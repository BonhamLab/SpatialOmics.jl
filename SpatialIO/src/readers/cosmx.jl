# SpatialIO/src/readers/cosmx.jl
# Reader for NanoString/Bruker CosMx SMI raw export data.
using Statistics: mean
#
# Raw export layout (slide-level gzipped CSVs):
#
#   <root>/
#   └── flatFiles/
#       └── <sample_name>/                                       ← sample dir
#           ├── <sample>_tx_file.csv.gz                         # transcripts (global coords)
#           ├── <sample>_fov_positions_file.csv.gz              # FOV origins
#           ├── <sample>-polygons.csv.gz                        # cell polygon vertices (global coords)
#           ├── <sample>_metadata_file.csv.gz                   # per-cell metadata
#           └── <sample>_exprMat_file.csv.gz                    # expression matrix
#   └── RawFiles/
#       └── <sample_name>/<run_timestamp>/CellStatsDir/
#           └── Morphology2D/                                   # per-FOV TIF files
#               └── <run>_C902_P99_N99_F<N:05d>.TIF
#
# Coordinate notes:
#   - x_global_px, y_global_px are already stitched in the raw-export CSVs.
#   - FOV coordinate relationship: x_global = x_fov_pos + x_local
#                                  y_global = y_fov_pos - y_local  (y-flip)
#   - TIF row 1 corresponds to y_local=0, i.e. y_global = y_fov_pos.
#   - FOV positions from fov_positions_file give the top-left corner of each FOV:
#       x_origin = x_global_px  (global x of TIF col 1)
#       y_origin = y_global_px  (global y of TIF row 1)

"""
    CosMxReader

Reader for NanoString/Bruker CosMx SMI raw export output.

Fields
------
- `lazy`             : Bool — reserved for future lazy transcript loading (default true)
- `sample`           : String or nothing — sample name (e.g. "mw_mus_p1_11"); auto-detected if nothing
- `morphology_dir`   : String or nothing — path to Morphology2D/ directory containing per-FOV TIFs
- `morphology_channel` : Int — which TIF channel to load (1-based; Morphology2D TIFs are 5-channel)
- `morphology_zarr`   : String or nothing — path where the morphology zarr store is written.
                        Defaults to a `morphology_cache.zarr` directory next to `morphology_dir`.
                        If the zarr already exists it is opened directly (no re-transcoding).

# Example
```julia
import SpatialIO as SIO
ds = SIO.load(SIO.CosMxReader(), "/data/my_experiment")
ds = SIO.load(SIO.CosMxReader(; morphology_dir="/data/RawFiles/Morphology2D"), "/data/my_experiment")
ds = SIO.load(SIO.CosMxReader(; morphology_dir="/data/RawFiles/Morphology2D",
                                morphology_zarr="/cache/morphology.zarr"), "/data/my_experiment")
```
"""
Base.@kwdef struct CosMxReader <: PlatformReader
    lazy::Bool = true
    sample::Union{Nothing, String} = nothing
    morphology_dir::Union{Nothing, String} = nothing
    morphology_channel::Int = 1
    morphology_zarr::Union{Nothing, String} = nothing
end

# ─── Path validation ──────────────────────────────────────────────────────────

function validate_path(reader::CosMxReader, path::String)::Bool
    isdir(path) || return false
    # Direct sample dir: contains *_tx_file.csv.gz
    any(endswith(f, "_tx_file.csv.gz") for f in readdir(path)) && return true
    # Parent with flatFiles/ subdir
    ff = joinpath(path, "flatFiles")
    isdir(ff) || return false
    # At least one sample subdir exists
    return any(isdir(joinpath(ff, d)) for d in readdir(ff))
end

# ─── Sample directory discovery ───────────────────────────────────────────────

"""
    _find_cosmx_sample(path, sample_name) -> sample_dir

Walk from `path` to the CosMx sample directory. Returns the directory that
directly contains the gzipped CSVs.

Accepts:
  - sample dir itself (contains `*_tx_file.csv.gz`)
  - parent dir containing `flatFiles/<sample>/`
"""
function _find_cosmx_sample(path::String, sample_name::Union{Nothing,String})
    # Already in a sample dir?
    if any(endswith(f, "_tx_file.csv.gz") for f in readdir(path))
        return path
    end
    ff = joinpath(path, "flatFiles")
    isdir(ff) || error("CosMxReader: cannot locate sample directory in $path")
    samples = sort(filter(d -> isdir(joinpath(ff, d)), readdir(ff)))
    isempty(samples) && error("CosMxReader: no sample subdirs in $ff")
    if sample_name !== nothing
        sample_name in samples ||
            error("CosMxReader: sample '$sample_name' not found in $ff; available: $samples")
        return joinpath(ff, sample_name)
    end
    length(samples) > 1 &&
        @warn "CosMxReader: multiple samples found in $ff, using '$(samples[1])'"
    return joinpath(ff, samples[1])
end

# ─── File locators ────────────────────────────────────────────────────────────

function _locate_cosmx_files(sample_dir::String)
    # Exclude macOS AppleDouble resource-fork files (._*) which are metadata-only
    files = filter(f -> !startswith(f, "._"), readdir(sample_dir))
    _find(suffix) = begin
        idx = findfirst(f -> endswith(f, suffix), files)
        isnothing(idx) && error("CosMxReader: no file ending '$suffix' in $sample_dir")
        joinpath(sample_dir, files[idx])
    end
    return (
        tx        = _find("_tx_file.csv.gz"),
        fov_pos   = _find("_fov_positions_file.csv.gz"),
        polygons  = _find("-polygons.csv.gz"),
        metadata  = _find("_metadata_file.csv.gz"),
        expr_mat  = _find("_exprMat_file.csv.gz"),
    )
end

# ─── FOV positions ────────────────────────────────────────────────────────────

function _read_fov_positions(path::String)
    df = CSV.read(path, DataFrame;
                  types = Dict("FOV" => Int32, "x_global_px" => Int32, "y_global_px" => Int32))
    return Dict{Int, NamedTuple{(:x, :y), Tuple{Int, Int}}}(
        row.FOV => (x = Int(row.x_global_px), y = Int(row.y_global_px))
        for row in eachrow(df)
    )
end

# ─── Transcripts ─────────────────────────────────────────────────────────────

function _read_transcripts(path::String)
    tx = CSV.read(path, DataFrame;
                  select = ["fov", "cell_ID", "cell", "x_global_px", "y_global_px", "z", "target", "CellComp"],
                  types  = Dict("fov" => Int32, "cell_ID" => Int32,
                                "x_global_px" => Float32, "y_global_px" => Float32,
                                "z" => Int8))
    coords = Matrix{Float32}(hcat(tx.x_global_px, tx.y_global_px))
    feats  = select(tx, Not([:x_global_px, :y_global_px]))
    return SpatialPoints(coords, feats, Dict{String,Any}("coord_cols" => ["x", "y"]))
end

# ─── Cell polygons ────────────────────────────────────────────────────────────

function _read_polygons(path::String)
    poly_df = CSV.read(path, DataFrame;
                       select = ["fov", "cellID", "cell", "x_global_px", "y_global_px"],
                       types  = Dict("fov" => Int32, "cellID" => Int32,
                                     "x_global_px" => Float64, "y_global_px" => Float64))

    grouped    = groupby(poly_df, :cell)
    n_cells    = length(grouped)
    geometries = Vector{Vector{UInt8}}(undef, n_cells)
    obs        = DataFrame(cell      = String[],
                           fov       = Int32[],
                           cell_ID   = Int32[],
                           x_centroid = Float64[],
                           y_centroid = Float64[])

    for (k, (key, grp)) in enumerate(pairs(grouped))
        xs = grp.x_global_px
        ys = grp.y_global_px
        ring = Point2f.(xs, ys)
        # Ensure ring is closed for WKB (first == last)
        if ring[1] != ring[end]
            push!(ring, ring[1])
        end
        poly          = GeometryBasics.Polygon(ring)
        geometries[k] = _encode_wkb_polygon(poly)
        push!(obs, (cell       = string(key.cell),
                    fov        = grp.fov[1],
                    cell_ID    = grp.cellID[1],
                    x_centroid = mean(xs),
                    y_centroid = mean(ys)))
    end

    return SpatialShapes(geometries, obs, Dict{String,Any}())
end

# ─── Expression matrix + metadata → SpatialTable ─────────────────────────────

function _read_expression(expr_path::String, meta_path::String)
    meta = CSV.read(meta_path, DataFrame)
    expr = CSV.read(expr_path, DataFrame;
                    types = Dict("fov" => Int32, "cell_ID" => Int32))

    # Join key: "fov:cell_ID"
    meta_key = string.(meta.fov, ":", meta.cell_ID)
    expr_key = string.(expr.fov, ":", expr.cell_ID)

    # Align expression rows to metadata order
    order = indexin(meta_key, expr_key)
    any(isnothing, order) &&
        error("CosMxReader: $(count(isnothing, order)) metadata cells missing from expression matrix")
    expr_sorted = expr[order, :]

    gene_names = names(expr)[3:end]   # skip fov, cell_ID
    X_dense    = Matrix{Float32}(expr_sorted[:, gene_names])
    X_sparse   = sparse(X_dense)
    obs = select(meta, Not(intersect(names(meta), ["_key"])))
    var = DataFrame(gene = gene_names)

    return SpatialTable(X_sparse, obs, var,
                        Dict{String,Any}("instance_key" => "cell",
                                         "region"       => "cell_boundaries"))
end

# ─── FOV grid → SpatialShapes ────────────────────────────────────────────────
#
# Stores each FOV as a rectangular polygon.  The features DataFrame carries
# pre-computed bbox columns so extent(ds.shapes["fovs"], i) is a fast lookup.

function _build_fov_shapes(
    fov_dict::Dict{Int, NamedTuple{(:x,:y),Tuple{Int,Int}}},
    fov_h::Int,
    fov_w::Int,
)::SpatialShapes
    fov_ids = sort(collect(keys(fov_dict)))
    geoms   = Vector{Any}(undef, length(fov_ids))
    obs     = DataFrame(
        fov_id     = Int32[],
        xmin       = Float64[],
        xmax       = Float64[],
        ymin       = Float64[],
        ymax       = Float64[],
        x_centroid = Float64[],
        y_centroid = Float64[],
    )
    for (k, fov_id) in enumerate(fov_ids)
        pos  = fov_dict[fov_id]
        xmin = Float64(pos.x)
        xmax = Float64(pos.x + fov_w - 1)
        ymax = Float64(pos.y)                # TIF row 1 = highest global y
        ymin = Float64(pos.y - fov_h + 1)
        ring = GeometryBasics.Point2f[
            (xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin),
        ]
        geoms[k] = GeometryBasics.Polygon(ring)
        push!(obs, (Int32(fov_id), xmin, xmax, ymin, ymax,
                    (xmin + xmax) / 2, (ymin + ymax) / 2))
    end
    return SpatialShapes(geoms, obs, Dict{String,Any}())
end

# ─── Morphology2D → zarr-backed SpatialImage ─────────────────────────────────
#
# Transcodes per-FOV TIF tiles into a single OME-NGFF zarr array on first load.
# Subsequent loads reopen the existing zarr directly — no re-transcoding.
#
# Array layout: (1, slide_height, slide_width) — OME-NGFF (c, y, x).
# Chunk layout: one chunk per FOV tile = (1, fov_height, fov_width).
# Y-flip: TIF row 1 (top of microscopy = highest global y) → last zarr row in
# the chunk, so zarr y index 0 = slide_ymin (bottom of slide).

function _write_morphology_zarr(
    zarr_path     ::String,
    morphology_dir::String,
    fov_dict      ::Dict{Int, NamedTuple{(:x,:y),Tuple{Int,Int}}},
    channel       ::Int,
)::SpatialImage
    # Fast path: zarr already built — open directly.
    if isdir(zarr_path) && isfile(joinpath(zarr_path, "zarr.json"))
        @info "CosMxReader: reusing cached morphology zarr" zarr_path
        return _read_sd_image(zarr_path)
    end

    @info "CosMxReader: transcoding Morphology2D TIFs to zarr" zarr_path channel

    # Discover TIF files: *_F<NNNNN>.TIF (case-insensitive)
    fov_paths = Dict{Int, String}()
    for f in readdir(morphology_dir)
        m = match(r"_F(\d+)\.TIF$"i, f)
        isnothing(m) && continue
        fov_paths[parse(Int, m[1])] = joinpath(morphology_dir, f)
    end
    isempty(fov_paths) &&
        error("CosMxReader: no TIF files found in $morphology_dir (pattern: *_F<N>.TIF)")

    # FOV dimensions from sample TIF (lazyio avoids loading pixels)
    sample_path = fov_paths[first(sort(collect(keys(fov_paths))))]
    sample_tif  = TiffImages.load(sample_path; lazyio = true)
    fov_h, fov_w = size(sample_tif, 1), size(sample_tif, 2)

    # Restrict to FOVs with both a TIF and a known grid position
    fov_ids   = sort(filter(f -> haskey(fov_dict, f), collect(keys(fov_paths))))
    isempty(fov_ids) &&
        error("CosMxReader: no overlap between TIF files and fov_positions — check morphology_dir")
    x_origins = Int[fov_dict[f].x for f in fov_ids]
    y_origins = Int[fov_dict[f].y for f in fov_ids]

    slide_xmin   = minimum(x_origins)
    slide_xmax   = maximum(x_origins) + fov_w - 1
    slide_ymin   = minimum(y_origins) - fov_h + 1
    slide_ymax   = maximum(y_origins)
    slide_width  = slide_xmax - slide_xmin + 1
    slide_height = slide_ymax - slide_ymin + 1

    # ── Write zarr store ──────────────────────────────────────────────────────

    mkpath(zarr_path)
    _zwrite_group(zarr_path, Dict{String,Any}(
        "x_range" => [slide_xmin, slide_xmax],
        "y_range" => [slide_ymin, slide_ymax],
    ))

    array_path  = joinpath(zarr_path, "0")
    chunk_shape = (1, fov_h, fov_w)
    mkpath(joinpath(array_path, "c"))
    Base.write(joinpath(array_path, "zarr.json"), JSON3.write(Dict(
        "zarr_format"  => 3,
        "node_type"    => "array",
        "shape"        => [1, slide_height, slide_width],
        "data_type"    => "float32",
        "chunk_grid"   => Dict("name" => "regular",
                               "configuration" => Dict("chunk_shape" => collect(chunk_shape))),
        "chunk_key_encoding" => Dict("name" => "default",
                                     "configuration" => Dict("separator" => "/")),
        "fill_value"   => 0,
        "codecs"       => [
            Dict("name" => "bytes", "configuration" => Dict("endian" => "little")),
            Dict("name" => "zstd",  "configuration" => Dict("level" => 0, "checksum" => false)),
        ],
        "attributes"   => Dict{String,Any}(),
        "storage_transformers" => [],
    )))

    for (fi, fov_id) in enumerate(fov_ids)
        xo = x_origins[fi]
        yo = y_origins[fi]
        chunk_x = div(xo - slide_xmin, fov_w)
        chunk_y = div(yo - fov_h + 1 - slide_ymin, fov_h)

        tif = TiffImages.load(fov_paths[fov_id]; lazyio = true)
        raw = ndims(tif) == 3 ? tif[:, :, channel] : tif[:, :]
        tile_flip = reverse(Float32.(raw); dims = 1)   # y-flip: TIF top → zarr bottom
        tile_3d   = reshape(tile_flip, 1, fov_h, fov_w)
        raw_bytes = reinterpret(UInt8, vec(permutedims(tile_3d, (3, 2, 1))))
        compressed = _zcompress(collect(raw_bytes))

        cpath = _chunk_path(array_path, (0, chunk_y, chunk_x))
        mkpath(dirname(cpath))
        Base.write(cpath, compressed)
    end

    return _read_sd_image(zarr_path)
end

# ─── Main read_data ───────────────────────────────────────────────────────────

function read_data(reader::CosMxReader, path::String)::SpatialDataset
    sample_dir = _find_cosmx_sample(path, reader.sample)
    sample_name = basename(sample_dir)
    @info "CosMxReader: loading sample" sample=sample_name sample_dir

    files    = _locate_cosmx_files(sample_dir)
    fov_dict = _read_fov_positions(files.fov_pos)

    @info "CosMxReader: reading transcripts"
    pts = _read_transcripts(files.tx)

    @info "CosMxReader: building cell polygons"
    shp = _read_polygons(files.polygons)

    @info "CosMxReader: reading expression matrix and cell metadata"
    tbl = _read_expression(files.expr_mat, files.metadata)

    ds = spatial_dataset(; metadata = Dict{String,Any}(
        "format"      => "CosMx",
        "sample"      => sample_name,
        "sample_dir"  => sample_dir,
        "fov_count"   => length(fov_dict),
    ))

    ds["transcripts"]    = pts
    ds["cell_boundaries"] = shp
    ds["expression"]     = tbl

    if reader.morphology_dir !== nothing
        zarr_path = something(reader.morphology_zarr,
                              abspath(joinpath(reader.morphology_dir, "..", "morphology_cache.zarr")))
        img = _write_morphology_zarr(zarr_path, reader.morphology_dir, fov_dict, reader.morphology_channel)
        ds["morphology"] = img
        # chunk shape is (c=1, y=fov_h, x=fov_w) — one chunk per FOV tile
        fov_h, fov_w = img.data.chunk_shape[2], img.data.chunk_shape[3]
        ds["fovs"] = _build_fov_shapes(fov_dict, fov_h, fov_w)
    end

    return ds
end
