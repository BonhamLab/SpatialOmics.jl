# SpatialIO/src/readers/cosmx.jl
# Reader for NanoString/Bruker CosMx SMI data in the DecodedFiles export layout.
#
# DecodedFiles layout (as of 2024-era exports):
#
#   <run_name>/
#   └── DecodedFiles/
#       └── <slide_name>/
#           └── <run_timestamp>/
#               ├── plex-<id>.txt                     # gene panel (DisplayName, CodeClass, ProbeID)
#               └── CellStatsDir/
#                   ├── CellComposite/
#                   │   └── CellComposite_F<N>.jpg    # composite images per FOV
#                   ├── CellOverlay/
#                   │   └── CellOverlay_F<N>.jpg
#                   ├── Morphology2D/                 # may be empty
#                   ├── Segmentation_<uuid>_<N>/      # alternative segmentation results (skip)
#                   └── FOV<N:05d>/                   # primary segmentation, one dir per FOV
#                       ├── CellBoundaries_F<N>.csv   # fov,cellID,x_local,y_local
#                       ├── CellBoundaries_F<N>_<cellID>.fz   # binary; not parsed
#                       ├── CellLabels_F<N>.tif       # cell label raster (Gray{N0f16}, 4256×4256)
#                       ├── CompartmentLabels_F<N>.tif
#                       ├── Run_<uuid>_<timestamp>_Cell_Stats_F<N>.csv
#                       │     columns: CellId, Area, AspectRatio, CenterX, CenterY,
#                       │              Width, Height, Mean-<marker>, Max-<marker>, ...
#                       └── Run_<uuid>_FOV<N:05d>__complete_code_cell_target_call_coord.csv
#                             columns: CellComp, CellId, ..., codeclass, fov,
#                                      seed_x, seed_y, target, x, y, z

# Coordinate notes:
#   - x, y in transcript CSV are FOV-local float pixel coordinates (0–4255.9).
#   - seed_x, seed_y are 10× sub-pixel coordinates relative to x, y.
#   - z = -1 indicates the projection plane; 0–7 are z-stack planes.
#   - No global positions file is present in this export format; global
#     stitching must be deferred or computed externally.
#   - FOV numbering is non-contiguous; always discover FOVs by scanning FOV* dirs.
#   - FOV image size is read from CellLabels TIF header (do not hardcode).

"""
    CosMxReader

Reader for NanoString/Bruker CosMx SMI output in the DecodedFiles export format.

Fields
------
- `lazy`        : Bool — defer loading large arrays until accessed (default true)
- `fov_subset`  : Union{Nothing, Vector{Int}} — restrict to specific FOV numbers,
                  or `nothing` to read all discovered FOVs (default nothing)
"""
Base.@kwdef struct CosMxReader <: PlatformReader
    lazy::Bool = true
    fov_subset::Union{Nothing, Vector{Int}} = nothing
end

"""
    load_cosmx(path; lazy=true, fov_subset=nothing) -> SpatialDataset

High-level convenience wrapper around `CosMxReader`.

`path` should point to the run directory containing `DecodedFiles/`, or any
directory in the hierarchy down to and including `CellStatsDir`.
"""
function load_cosmx(
    path::String;
    lazy::Bool = true,
    fov_subset::Union{Nothing, Vector{Int}} = nothing,
)
    return load_spatial_data(path, CosMxReader; lazy = lazy, fov_subset = fov_subset)
end

"""
    validate_path(reader::CosMxReader, path::String) -> Bool

Return `true` if `path` looks like a CosMx DecodedFiles export.
"""
function validate_path(reader::CosMxReader, path::String)::Bool
    isdir(path) || return false
    basename(path) == "CellStatsDir" && return true
    basename(path) == "DecodedFiles" && return true
    return isdir(joinpath(path, "DecodedFiles"))
end

# ─── Directory walking ────────────────────────────────────────────────────────

"""
    find_cell_stats_dir(root::String) -> (cell_stats_dir, run_dir)

Walk from `root` to `CellStatsDir` and return both the CellStatsDir path and
the run-timestamp directory (parent of CellStatsDir, which contains plex-*.txt).

Accepts any of:
  - the run name dir containing DecodedFiles/
  - the DecodedFiles/ dir itself
  - the slide dir
  - the run-timestamp dir
  - CellStatsDir itself
"""
function find_cell_stats_dir(root::String)
    isdir(root) || error("find_cell_stats_dir: not a directory: $root")

    # Already at CellStatsDir
    if basename(root) == "CellStatsDir"
        return root, dirname(root)
    end

    # Descend through DecodedFiles → slide → run → CellStatsDir
    df_root = if basename(root) == "DecodedFiles"
        root
    elseif isdir(joinpath(root, "DecodedFiles"))
        joinpath(root, "DecodedFiles")
    else
        # Maybe root IS slide or run dir — try looking for CellStatsDir directly
        csd = joinpath(root, "CellStatsDir")
        isdir(csd) && return csd, root
        # One level up: root is slide dir
        runs = _sorted_subdirs(root)
        for r in reverse(runs)
            csd = joinpath(root, r, "CellStatsDir")
            isdir(csd) && return csd, joinpath(root, r)
        end
        error("find_cell_stats_dir: cannot locate CellStatsDir from $root")
    end

    slides = _sorted_subdirs(df_root)
    isempty(slides) && error("find_cell_stats_dir: no slide dirs in $df_root")
    length(slides) > 1 &&
        @warn "find_cell_stats_dir: multiple slides found, using $(slides[1])"
    slide_path = joinpath(df_root, slides[1])

    runs = _sorted_subdirs(slide_path)
    isempty(runs) && error("find_cell_stats_dir: no run dirs in $slide_path")
    length(runs) > 1 &&
        @warn "find_cell_stats_dir: multiple runs found, using most recent: $(last(runs))"
    run_path = joinpath(slide_path, last(runs))

    csd = joinpath(run_path, "CellStatsDir")
    isdir(csd) || error("find_cell_stats_dir: CellStatsDir not found in $run_path")
    return csd, run_path
end

_sorted_subdirs(dir) = sort(filter(d -> isdir(joinpath(dir, d)), readdir(dir)))

"""
    discover_fovs(cell_stats_dir::String) -> Vector{Int}

Scan `cell_stats_dir` for FOV<N> subdirectories and return sorted FOV numbers.
Non-contiguous FOV numbering is expected and handled correctly.
"""
function discover_fovs(cell_stats_dir::String)::Vector{Int}
    skip = Set(["CellComposite", "CellOverlay", "Morphology2D"])
    fov_ids = Int[]
    for d in readdir(cell_stats_dir)
        d in skip && continue
        startswith(d, "Segmentation_") && continue
        m = match(r"^FOV(\d+)$", d)
        isnothing(m) && continue
        push!(fov_ids, parse(Int, m[1]))
    end
    return sort(fov_ids)
end

# ─── Per-FOV file locators ────────────────────────────────────────────────────

function _fov_dir(cell_stats_dir, fov_id)
    joinpath(cell_stats_dir, @sprintf("FOV%05d", fov_id))
end

function _locate_transcript_csv(fov_dir, fov_id)
    pattern = Regex("^Run_.*_FOV$(lpad(fov_id, 5, '0'))__complete_code_cell_target_call_coord\\.csv\$")
    f = _find_file(fov_dir, pattern)
    if isnothing(f)
        @warn "No transcript CSV found for FOV $fov_id in $fov_dir — skipping transcripts for this FOV"
    end
    return f
end

function _locate_cell_stats_csv(fov_dir, fov_id)
    pattern = Regex("^Run_.*_Cell_Stats_F$(lpad(fov_id, 5, '0'))\\.csv\$")
    f = _find_file(fov_dir, pattern)
    isnothing(f) && error("Cell stats CSV not found in $fov_dir")
    return f
end

function _find_file(dir, pattern::Regex)
    for f in readdir(dir)
        occursin(pattern, f) && return joinpath(dir, f)
    end
    return nothing
end

# ─── Plex / gene panel ────────────────────────────────────────────────────────

function _find_plex_file(run_dir)
    for f in readdir(run_dir)
        startswith(f, "plex-") && endswith(f, ".txt") && return joinpath(run_dir, f)
    end
    return nothing
end

function _read_plex(path)
    isnothing(path) && return DataFrame(DisplayName=String[], CodeClass=String[], ProbeID=String[])
    return CSV.read(path, DataFrame)
end

# ─── TIF dimension reader ─────────────────────────────────────────────────────

"""Read only the (height, width) from a TIFF header without loading pixel data."""
function _tif_size(path)
    img = TiffImages.load(path; lazyio=true)
    return size(img, 1), size(img, 2)   # drop trailing singleton channel dim
end

# ─── Main read_data ───────────────────────────────────────────────────────────

function read_data(reader::CosMxReader, path::String)::SpatialDataset
    cell_stats_dir, run_dir = find_cell_stats_dir(path)

    all_fov_ids = discover_fovs(cell_stats_dir)
    fov_ids = if isnothing(reader.fov_subset)
        all_fov_ids
    else
        ids = intersect(reader.fov_subset, all_fov_ids)
        isempty(ids) && error("fov_subset $(reader.fov_subset) has no overlap with discovered FOVs $all_fov_ids")
        sort(ids)
    end

    plex_df = _read_plex(_find_plex_file(run_dir))

    # Accumulate per-FOV data
    tx_frames    = Vector{DataFrame}(undef, length(fov_ids))
    stats_frames = Vector{DataFrame}(undef, length(fov_ids))
    bounds_frames = Vector{DataFrame}(undef, length(fov_ids))
    fov_sizes    = Dict{Int, Tuple{Int,Int}}()   # fov_id => (h, w)

    for (i, fov_id) in enumerate(fov_ids)
        fdir = _fov_dir(cell_stats_dir, fov_id)

        tx_path = _locate_transcript_csv(fdir, fov_id)
        tx = if isnothing(tx_path)
            DataFrame(CellComp=String[], CellId=Int32[], codeclass=String[],
                      fov=Int32[], target=String[],
                      x=Float32[], y=Float32[], z=Int8[])
        else
            CSV.read(tx_path, DataFrame;
                     select = ["CellComp", "CellId", "codeclass", "fov",
                               "target", "x", "y", "z"],
                     types  = Dict("fov" => Int32, "CellId" => Int32,
                                   "x" => Float32, "y" => Float32,
                                   "z" => Int8))
        end
        tx_frames[i] = tx

        st = CSV.read(_locate_cell_stats_csv(fdir, fov_id), DataFrame)
        st[!, :fov] .= Int32(fov_id)
        stats_frames[i] = st

        bc = CSV.read(joinpath(fdir, "CellBoundaries_F$(lpad(fov_id,5,'0')).csv"), DataFrame)
        bounds_frames[i] = bc

        tif_path = joinpath(fdir, "CellLabels_F$(lpad(fov_id,5,'0')).tif")
        fov_sizes[fov_id] = _tif_size(tif_path)
    end

    all_tx    = vcat(tx_frames...)
    all_stats = vcat(stats_frames...)
    all_bounds = vcat(bounds_frames...)

    ds = spatial_dataset(; metadata = Dict{String,Any}(
        "format"    => "CosMx-DecodedFiles",
        "run_dir"   => run_dir,
        "fov_ids"   => fov_ids,
        "fov_sizes" => fov_sizes,
    ))

    # ── SpatialPoints — transcripts ───────────────────────────────────────────
    coords = Matrix{Float32}(hcat(all_tx.x, all_tx.y))   # N×2
    feat   = select(all_tx, Not([:x, :y]))
    add_points!(ds, "transcripts", SpatialPoints(coords, feat, Dict{String,Any}()))

    # ── SpatialTable — expression counts ─────────────────────────────────────
    # Build a global cell key "fov_cellid" to uniquify cells across FOVs
    endo_tx = filter(r -> r.codeclass == "Endogenous" && r.CellId > 0, all_tx)

    if !isempty(endo_tx) && !isempty(plex_df)
        expr_tbl = _build_expression_table(endo_tx, all_stats, plex_df)
        add_table!(ds, "table", expr_tbl)
    end

    # ── SpatialShapes — cell boundaries ──────────────────────────────────────
    if !isempty(all_bounds)
        shp = _build_shapes(all_bounds)
        add_shapes!(ds, "cell_boundaries", shp)
    end

    # ── SpatialImage — composite JPEGs (skip if lazy) ────────────────────────
    if !reader.lazy
        composite_dir = joinpath(cell_stats_dir, "CellComposite")
        if isdir(composite_dir)
            for fov_id in fov_ids
                jpg = joinpath(composite_dir, "CellComposite_F$(lpad(fov_id,5,'0')).jpg")
                isfile(jpg) || continue
                img = _load_composite_jpg(jpg)
                add_image!(ds, "composite_F$(lpad(fov_id,5,'0'))", img)
            end
        end
    end

    return ds
end

# ─── Expression count matrix ──────────────────────────────────────────────────

function _build_expression_table(
    endo_tx::DataFrame,
    stats::DataFrame,
    plex::DataFrame,
)::SpatialTable
    # Global cell index: "fov_cellid" string key → integer row index
    cell_keys = sort(unique(string.(endo_tx.fov) .* "_" .* string.(endo_tx.CellId)))
    cell_idx  = Dict(k => i for (i, k) in enumerate(cell_keys))
    n_cells   = length(cell_keys)

    # Gene index from plex (Endogenous only, preserving plex order)
    gene_names = plex[plex.CodeClass .== "Endogenous", :DisplayName]
    gene_idx   = Dict(g => i for (i, g) in enumerate(gene_names))
    n_genes    = length(gene_names)

    # Build COO lists
    I_vals = Int32[]
    J_vals = Int32[]
    V_vals = Float32[]

    for gdf in groupby(endo_tx, [:fov, :CellId, :target])
        key = string(gdf.fov[1]) * "_" * string(gdf.CellId[1])
        ci  = get(cell_idx, key, nothing)
        gi  = get(gene_idx, gdf.target[1], nothing)
        (isnothing(ci) || isnothing(gi)) && continue
        push!(I_vals, ci)
        push!(J_vals, gi)
        push!(V_vals, Float32(nrow(gdf)))
    end

    X = sparse(I_vals, J_vals, V_vals, n_cells, n_genes)

    # obs DataFrame — join cell metadata from stats
    obs = DataFrame(cell_key = cell_keys)
    obs[!, :fov]    = Int32.(parse.(Int, first.(split.(cell_keys, "_"))))
    obs[!, :CellId] = Int32.(parse.(Int, last.(split.(cell_keys, "_"))))
    stats_sub = select(stats, [:fov, :CellId, :Area, :CenterX, :CenterY])
    leftjoin!(obs, stats_sub; on = [:fov, :CellId])
    # Fill missing stats (cells present in transcripts but absent from cell stats)
    for col in (:Area, :CenterX, :CenterY)
        obs[!, col] = Int32.(coalesce.(obs[!, col], Int32(0)))
    end

    var = select(plex[plex.CodeClass .== "Endogenous", :], [:DisplayName, :CodeClass])
    rename!(var, :DisplayName => :gene)

    return SpatialTable(X, obs, var, Dict{String,Any}())
end

# ─── Cell boundary shapes ─────────────────────────────────────────────────────

function _build_shapes(bounds::DataFrame)::SpatialShapes
    # Group boundary vertices per (fov, cellID) → store as raw vertex matrix
    # WKB encoding is deferred; geometries are Vector{Matrix{Int}} for now
    geoms = Vector{Any}()
    keys_df = unique(select(bounds, [:fov, :cellID]))
    for row in eachrow(keys_df)
        sub = bounds[(bounds.fov .== row.fov) .& (bounds.cellID .== row.cellID), :]
        push!(geoms, Matrix{Int}(hcat(sub.x_local, sub.y_local)))
    end
    feat = copy(keys_df)
    return SpatialShapes(geoms, feat, Dict{String,Any}())
end

# ─── Composite image loader ───────────────────────────────────────────────────

function _load_composite_jpg(path::String)::SpatialImage
    # FileIO is not a SpatialIO dep; use a basic JPEG read via ImageIO if available,
    # otherwise store the path in metadata for later loading.
    # TODO: add ImageIO as an optional dep and load lazily
    return SpatialImage(
        zeros(UInt8, 0, 0, 0),          # placeholder — not loaded
        (;),
        Dict{String,Any}("jpeg_path" => path),
    )
end
