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
- `lazy`           : Bool — reserved for future lazy transcript loading (default true)
- `sample`         : String or nothing — sample name (e.g. "mw_mus_p1_11"); auto-detected if nothing
- `morphology_dir` : String or nothing — path to Morphology2D/ directory containing per-FOV TIFs.
                     All channels are transcoded; channel names are read from TIF metadata and stored
                     in `img.metadata["channel_names"]`.
- `morphology_zarr` : String or nothing — path where the morphology zarr store is written.
                      Defaults to `morphology_cache.zarr` next to `morphology_dir`.
                      Reopened on subsequent loads; rebuilt automatically if the channel count changes.

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
    morphology_zarr::Union{Nothing, String} = nothing
end

function Base.show(io::IO, r::CosMxReader)
    parts = String[]
    r.sample         !== nothing && push!(parts, "sample=$(repr(r.sample))")
    r.morphology_dir !== nothing && push!(parts, "morphology_dir=$(repr(r.morphology_dir))")
    r.morphology_zarr !== nothing && push!(parts, "morphology_zarr=$(repr(r.morphology_zarr))")
    print(io, "CosMxReader(", join(parts, ", "), ")")
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
    # Normalize to canonical schema (see platform_reader.jl)
    rename!(tx, :target => :gene, :cell_ID => :cell_id, :CellComp => :cell_compartment)
    return SpatialPoints(tx; x_col=:x_global_px, y_col=:y_global_px, label_col=:gene)
end

# ─── Cell polygons ────────────────────────────────────────────────────────────

function _read_polygons(path::String)
    poly_df = CSV.read(path, DataFrame;
                       select = ["fov", "cellID", "cell", "x_global_px", "y_global_px"],
                       types  = Dict("fov" => Int32, "cellID" => Int32,
                                     "x_global_px" => Float64, "y_global_px" => Float64))

    grouped = groupby(poly_df, :cell)
    shapes  = [begin
        xs   = grp.x_global_px
        ys   = grp.y_global_px
        ring = Point2f.(xs, ys)
        ring[1] != ring[end] && push!(ring, ring[1])
        poly = GeometryBasics.Polygon(ring)
        SpatialShape(poly, string(key.cell), (cell = string(key.cell), cell_ID = grp.cellID[1]))
    end for (key, grp) in pairs(grouped)]

    return SpatialShapes(shapes, Dict{String,Any}())
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
    shapes  = [begin
        pos  = fov_dict[fov_id]
        xmin = Float64(pos.x)
        xmax = Float64(pos.x + fov_w - 1)
        ymax = Float64(pos.y)
        ymin = Float64(pos.y - fov_h + 1)
        ring = GeometryBasics.Point2f[
            (xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin),
        ]
        SpatialShape(GeometryBasics.Polygon(ring), string(fov_id), (fov_id = Int32(fov_id),))
    end for fov_id in fov_ids]
    return SpatialShapes(shapes, Dict{String,Any}())
end

# ─── Morphology2D → zarr-backed SpatialImage ─────────────────────────────────
#
# Transcodes per-FOV TIF tiles into a single OME-NGFF zarr array on first load.
# ALL channels are written in one pass; subsequent loads reopen the zarr directly.
#
# Array layout: (n_channels, slide_height, slide_width) — OME-NGFF (c, y, x).
# Chunk shape: (n_channels, fov_h, fov_w).  FOV origins are NOT assumed to be
# on a multiple-of-fov_w/fov_h grid — each FOV is placed at its exact pixel
# offset.  A partial FOV tile may therefore span up to four chunk boundaries.
# The implementation builds a full slide buffer in RAM, fills tiles at their
# exact positions, then writes chunks in a single pass.
# Y-flip: TIF row 1 (top of microscopy = highest global y) → last slide row
# in that FOV's pixel range, so slide row 0 = slide_ymin (bottom of slide).
#
# Channel names are parsed from the IMAGEDESCRIPTION JSON in the first TIF IFD
# and stored in img.metadata["channel_names"].  The cache is invalidated and
# rebuilt when n_channels changes or the "pixel_exact" flag is absent (caches
# written by the old chunk-snapping code lack this flag).

"""
    _parse_channel_names(tif_path) -> Vector{String}

Extract biological target names from the CosMx TIF IMAGEDESCRIPTION JSON.
Returns e.g. `["PanCK", "CD68", "CD298/B2M", "CD45", "DAPI"]` ordered by
`ChannelOrder` (BGYRU).  Falls back to `["channel_1", …]` on any parse error.
"""
function _parse_channel_names(tif_path::String)::Vector{String}
    tif = TiffImages.load(tif_path; lazyio = true)
    n   = ndims(tif) == 3 ? size(tif, 3) : 1
    fallback = ["channel_$i" for i in 1:n]

    desc = nothing
    for (tag, val) in first(TiffImages.ifds(tif))
        if tag == 270      # IMAGEDESCRIPTION
            desc = val[1].data
            break
        end
    end
    desc === nothing && return fallback

    meta = try JSON.parse(desc, Dict) catch; return fallback end
    channel_order = get(meta, "ChannelOrder", nothing)
    reagents      = get(get(meta, "MorphologyKit", Dict()), "MorphologyReagents", nothing)
    (channel_order === nothing || reagents === nothing) && return fallback

    id_to_target = Dict{String,String}()
    for r in reagents
        flu = get(r, "Fluorophore", nothing)
        flu === nothing && continue
        cid = get(flu, "ChannelId", nothing)
        tgt = get(r, "BiologicalTarget", nothing)
        (cid !== nothing && tgt !== nothing) && (id_to_target[string(cid)] = string(tgt))
    end

    return [get(id_to_target, string(c), "channel_$i")
            for (i, c) in enumerate(channel_order)]
end

function _write_morphology_zarr(
    zarr_path     ::String,
    morphology_dir::String,
    fov_dict      ::Dict{Int, NamedTuple{(:x,:y),Tuple{Int,Int}}},
)::SpatialImage
    # Discover TIF files: *_F<NNNNN>.TIF (case-insensitive)
    fov_paths = Dict{Int, String}()
    for f in readdir(morphology_dir)
        m = match(r"_F(\d+)\.TIF$"i, f)
        isnothing(m) && continue
        fov_paths[parse(Int, m[1])] = joinpath(morphology_dir, f)
    end
    isempty(fov_paths) &&
        error("CosMxReader: no TIF files found in $morphology_dir (pattern: *_F<N>.TIF)")

    # Peek at one TIF to get dimensions and channel count (lazyio = header only)
    sample_path  = fov_paths[first(sort(collect(keys(fov_paths))))]
    sample_tif   = TiffImages.load(sample_path; lazyio = true)
    fov_h, fov_w = size(sample_tif, 1), size(sample_tif, 2)
    n_channels   = ndims(sample_tif) == 3 ? size(sample_tif, 3) : 1
    channel_names = _parse_channel_names(sample_path)

    # Fast path: reuse existing zarr if channel count matches AND was written
    # with exact pixel placement (pixel_exact flag in root zarr.json).
    array_meta_path = joinpath(zarr_path, "0", "zarr.json")
    root_meta_path  = joinpath(zarr_path, "zarr.json")
    if isdir(zarr_path) && isfile(array_meta_path)
        existing     = JSON.parse(Base.read(array_meta_path, String), Dict)
        root_attrs   = if isfile(root_meta_path)
            get(JSON.parse(Base.read(root_meta_path, String), Dict), "attributes", Dict())
        else
            Dict()
        end
        pixel_exact  = get(root_attrs, "pixel_exact", false)
        if existing["shape"][1] == n_channels && pixel_exact && isdir(joinpath(zarr_path, "1"))
            @info "CosMxReader: reusing cached morphology zarr" zarr_path
            return _read_sd_image(zarr_path)
        elseif existing["shape"][1] != n_channels
            @info "CosMxReader: channel count changed — rebuilding zarr" zarr_path old=existing["shape"][1] new=n_channels
            rm(zarr_path; recursive=true)
        else
            @info "CosMxReader: rebuilding zarr with exact pixel placement" zarr_path
            rm(zarr_path; recursive=true)
        end
    end

    @info "CosMxReader: transcoding Morphology2D TIFs to zarr" zarr_path n_channels channel_names

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
        "x_range"       => [slide_xmin, slide_xmax],
        "y_range"       => [slide_ymin, slide_ymax],
        "channel_names" => channel_names,
        "pixel_exact"   => true,
    ))

    array_path  = joinpath(zarr_path, "0")
    mkpath(joinpath(array_path, "c"))
    Base.write(joinpath(array_path, "zarr.json"), JSON.json(Dict(
        "zarr_format"  => 3,
        "node_type"    => "array",
        "shape"        => [n_channels, slide_height, slide_width],
        "data_type"    => "float32",
        "chunk_grid"   => Dict("name" => "regular",
                               "configuration" => Dict("chunk_shape" => [1, fov_h, fov_w])),
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

    # Level-1: coarse pyramid level for zoom-out display via ImagePyramidSampler.
    # Without a pyramid the heatmap path falls back to strided zarr reads over the
    # full ~100k×90k slide, which either OOMs or returns zeros.
    ds         = 32
    l1_h       = cld(fov_h, ds)
    l1_w       = cld(fov_w, ds)
    l1_slide_h = cld(slide_height, ds)
    l1_slide_w = cld(slide_width,  ds)
    level1_path = joinpath(zarr_path, "1")
    mkpath(joinpath(level1_path, "c"))
    Base.write(joinpath(level1_path, "zarr.json"), JSON.json(Dict(
        "zarr_format"  => 3,
        "node_type"    => "array",
        "shape"        => [n_channels, l1_slide_h, l1_slide_w],
        "data_type"    => "float32",
        "chunk_grid"   => Dict("name" => "regular",
                               "configuration" => Dict("chunk_shape" => [1, l1_h, l1_w])),
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

    # ── Precompute FOV pixel ranges (0-based slide coordinates) ──────────────
    n_fovs  = length(fov_ids)
    fov_r1s = [y_origins[fi] - fov_h + 1 - slide_ymin for fi in 1:n_fovs]
    fov_r2s = [y_origins[fi]             - slide_ymin for fi in 1:n_fovs]
    fov_c1s = [x_origins[fi]             - slide_xmin for fi in 1:n_fovs]
    fov_c2s = [x_origins[fi] + fov_w - 1 - slide_xmin for fi in 1:n_fovs]

    n_cy = cld(slide_height, fov_h)
    n_cx = cld(slide_width,  fov_w)

    # ── Write one (channel, cy, cx) block at a time ───────────────────────────
    # Chunk shape (1, fov_h, fov_w): one zarr chunk per channel per spatial tile.
    # Peak memory per iteration: one 2-D float32 tile (~69 MB) + one TIF channel
    # materialised (~72 MB) — vs the previous (n_channels × fov_h × fov_w) buffer
    # that consumed ~345 MB.  The same TIF is loaded n_channels times per overlapping
    # FOV; subsequent loads hit the OS page cache so I/O cost is low.

    for cy in 0:(n_cy - 1), cx in 0:(n_cx - 1)
        # Chunk pixel range in slide (0-based, inclusive)
        cr1 = cy * fov_h;  cr2 = min(cr1 + fov_h - 1, slide_height - 1)
        cc1 = cx * fov_w;  cc2 = min(cc1 + fov_w - 1, slide_width  - 1)

        # Collect overlapping FOV indices
        overlapping = Int[]
        for fi in 1:n_fovs
            (fov_r1s[fi] <= cr2 && fov_r2s[fi] >= cr1 &&
             fov_c1s[fi] <= cc2 && fov_c2s[fi] >= cc1) && push!(overlapping, fi)
        end
        isempty(overlapping) && continue   # no file written; zarr fill_value = 0

        # Precompute src/dst index ranges once (independent of channel)
        fov_coords = [
            let ir1 = max(fov_r1s[fi], cr1), ir2 = min(fov_r2s[fi], cr2),
                ic1 = max(fov_c1s[fi], cc1), ic2 = min(fov_c2s[fi], cc2)
                (dst_r1 = ir1-cr1+1,          dst_r2 = ir2-cr1+1,
                 dst_c1 = ic1-cc1+1,          dst_c2 = ic2-cc1+1,
                 src_r1 = ir1-fov_r1s[fi]+1,  src_r2 = ir2-fov_r1s[fi]+1,
                 src_c1 = ic1-fov_c1s[fi]+1,  src_c2 = ic2-fov_c1s[fi]+1)
            end
            for fi in overlapping]

        for c in 1:n_channels
            chunk_c = zeros(Float32, fov_h, fov_w)

            for (fi, coord) in zip(overlapping, fov_coords)
                tif     = TiffImages.load(fov_paths[fov_ids[fi]]; lazyio = true)
                raw     = ndims(tif) == 3 ? tif[:, :, c] : tif[:, :]
                flipped = reverse(Float32.(raw); dims = 1)   # y-flip
                chunk_c[coord.dst_r1:coord.dst_r2, coord.dst_c1:coord.dst_c2] .=
                    flipped[coord.src_r1:coord.src_r2, coord.src_c1:coord.src_c2]
            end

            # Write level-0 chunk (channel c-1, spatial tile cy,cx)
            buf0       = reshape(chunk_c, 1, fov_h, fov_w)
            raw_bytes0 = reinterpret(UInt8, vec(permutedims(buf0, (3, 2, 1))))
            cpath0     = _chunk_path(array_path, (c-1, cy, cx))
            mkpath(dirname(cpath0))
            Base.write(cpath0, _zcompress(collect(raw_bytes0)))

            # Write level-1 chunk (ds× spatial subsampling)
            l1_raw = chunk_c[1:ds:end, 1:ds:end]
            rh, rw = size(l1_raw)
            buf1   = zeros(Float32, 1, l1_h, l1_w)
            buf1[1, 1:rh, 1:rw] .= l1_raw
            raw_bytes1 = reinterpret(UInt8, vec(permutedims(buf1, (3, 2, 1))))
            cpath1     = _chunk_path(level1_path, (c-1, cy, cx))
            mkpath(dirname(cpath1))
            Base.write(cpath1, _zcompress(collect(raw_bytes1)))
        end
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
        img = _write_morphology_zarr(zarr_path, reader.morphology_dir, fov_dict)
        ds["morphology"] = img
        # chunk shape is (n_channels, fov_h, fov_w) — one chunk per FOV tile
        fov_h, fov_w = img.data.chunk_shape[2], img.data.chunk_shape[3]
        ds["fovs"] = _build_fov_shapes(fov_dict, fov_h, fov_w)
    end

    return ds
end
