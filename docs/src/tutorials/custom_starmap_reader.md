# Build a custom STARmap reader

```@setup starmap-figure
using SpatialOmics
using CairoMakie
using Markdown
CairoMakie.activate!(type="svg")
set_theme!(Theme(
    fontsize=15,
    Figure=(; backgroundcolor=:white),
    Axis=(; xgridvisible=false, ygridvisible=false),
))
```

STARmap export layouts are still evolving. This tutorial shows how to adapt one observed
STARmap PLUS layout without making that layout part of SpatialOmics' stable API. The same
pattern applies to an internal assay or a new vendor format: use a small format token,
validate assumptions at the boundary, and return ordinary SpatialOmics elements.

The example layout contains:

- `Dapi/configurations.registered.txt`, with `(x, y, z)` offsets for each tile;
- `Dapi/*_Tile_NNN_405.tif`, with one DAPI z-stack per tile;
- `Decoded_spots/Raw/fused_goodSpots.csv`, in canvas-relative coordinates;
- a filtered fused transcript CSV; and
- optional per-tile transcript CSVs in tile-local coordinates.

Check these names and coordinate conventions against each new export. A reader should fail
clearly when its assumptions do not hold.

## Define a format token

Format tokens keep dispatch at the I/O boundary and avoid technology checks throughout an
analysis workflow.

```julia
using SpatialOmics
using CSV, DataFrames
import Images: load

struct STARmap end

function parse_registration(path::AbstractString)
    entries = Pair{String,NTuple{3,Float64}}[]
    for line in eachline(path)
        match_result = match(r"^(Tile_\d+\.tif);\s*;\s*\(([^)]+)\)", line)
        isnothing(match_result) && continue
        values = parse.(Float64, split(match_result[2], ","))
        length(values) == 3 || error("Expected three offsets in: $line")
        push!(entries, match_result[1] => (values[1], values[2], values[3]))
    end
    isempty(entries) && error("No tile registrations found in $path")
    entries
end
```

Natural errors from `CSV.read`, `load`, and constructors are useful here. Avoid a broad
`try`/`catch` that turns a malformed export into a partially populated dataset.

## Preserve raw columns while defining a working coordinate system

In this export, fused transcript coordinates start at the stitched canvas origin while
registered tile offsets may be negative. Working coordinates are shifted into the
registration frame, but the original columns are retained so an exporter can reconstruct
the vendor table.

```julia
function transcript_points(table, x_offset, y_offset)
    canvas_x = copy(table.x)
    canvas_y = copy(table.y)
    table.x .+= x_offset
    table.y .+= y_offset
    SpatialPoints(
        table;
        gene=:gene,
        features=(z=table.z, canvas_x, canvas_y),
        coord_system="global_px",
    )
end
```

Keeping raw columns is different from keeping a reference to the raw files. For lossless
recovery, record both.

## Construct elements and acquisition sources

The fused CSV does not identify which tile produced a transcript. Do not infer a unique
origin from position: a transcript can lie in two overlapping tile footprints. Footprints
can still carry exact source provenance, and each source can retain paths and registration
values needed to reopen the raw data.

```julia
function Base.read(::STARmap, root::AbstractString)
    registration_file = joinpath(root, "Dapi", "configurations.registered.txt")
    registration = parse_registration(registration_file)

    dapi_files = filter(
        path -> endswith(basename(path), "_405.tif"),
        readdir(joinpath(root, "Dapi"); join=true),
    )
    isempty(dapi_files) && error("No DAPI TIFFs found under $root")
    tile_height, tile_width = size(load(first(dapi_files)))[1:2]

    x_offset = floor(Int, minimum(value[1] for (_, value) in registration))
    y_offset = floor(Int, minimum(value[2] for (_, value) in registration))
    raw_path = joinpath(root, "Decoded_spots", "Raw", "fused_goodSpots.csv")
    filtered_path = joinpath(
        root,
        "Decoded_spots",
        "Cr2_Cxcl13_clusters_removed_fused_goodSpots.csv",
    )

    dataset = SpatialDataset(metadata=Dict(
        "technology" => "STARmap PLUS",
        "source_root" => abspath(root),
        "raw_transcripts" => relpath(raw_path, root),
        "filtered_transcripts" => relpath(filtered_path, root),
    ))
    push!(dataset, CoordinateSystem("global_px"; units=("px", "px")))
    dataset["transcripts_raw"] = transcript_points(
        CSV.read(raw_path, DataFrame), x_offset, y_offset,
    )
    dataset["transcripts_filtered"] = transcript_points(
        CSV.read(filtered_path, DataFrame), x_offset, y_offset,
    )

    tile_names = [splitext(name)[1] for (name, _) in registration]
    footprints = [
        Polygon(Point2f[
            (x, y), (x + tile_width, y),
            (x + tile_width, y + tile_height), (x, y + tile_height), (x, y),
        ])
        for (_, (x, y, _)) in registration
    ]
    dataset["tile_footprints"] = SpatialShapes(
        footprints;
        instance_id=Int32.(eachindex(footprints)),
        origins=tile_names,
        coord_system="global_px",
    )

    for (index, ((registered_name, offset), tile_name)) in
            enumerate(zip(registration, tile_names))
        tile_number = match(r"Tile_(\d+)", registered_name)[1]
        dapi_path = only(filter(
            path -> occursin("Tile_$(tile_number)_", basename(path)), dapi_files,
        ))
        push!(dataset, AcquisitionSource(
            tile_name;
            region="tile_footprints",
            instance_id=index,
            attributes=Dict(
                "technology" => "STARmap PLUS",
                "registration_name" => registered_name,
                "offset_xyz_px" => collect(offset),
                "dapi_path" => relpath(dapi_path, root),
                "raw_spots_path" => joinpath(
                    "Decoded_spots", "Raw", "$(tile_name)_goodSpots.csv",
                ),
            ),
        ))
    end
    dataset
end
```

Inspect the result through public accessors:

```julia
dataset = read(STARmap(), "/path/to/export")
sources(dataset)
source_attributes(dataset, first(sources(dataset)))
features(points(dataset, "transcripts_raw"), :canvas_x)
```

A source view of `tile_footprints` is provenance-exact. A source view of the fused
transcript element emits a warning and uses the source footprint because the CSV did not
provide transcript origins. Reading the per-tile CSVs is the right extension when exact
tile-level transcript provenance is required.

The geometry explains why position cannot recover a unique acquisition source:
transcripts in the overlap are compatible with both registered tile footprints.

```@eval starmap-figure
tile_1 = Polygon(Point2f[(0, 0), (8, 0), (8, 7), (0, 7), (0, 0)])
tile_2 = Polygon(Point2f[(5, 2), (13, 2), (13, 9), (5, 9), (5, 2)])
fused_transcripts = Point2f[(2, 2), (6, 3), (7, 5), (10, 7), (12, 4)]
figure = Figure(size=(720, 430))
axis = Axis(figure[1, 1]; aspect=DataAspect(), title="Registered STARmap tiles",
            xlabel="global x (px)", ylabel="global y (px)")
poly!(axis, [tile_1]; color=(:dodgerblue, 0.2), strokecolor=:dodgerblue3,
      strokewidth=2, label="Tile 1 footprint")
poly!(axis, [tile_2]; color=(:darkorange, 0.2), strokecolor=:darkorange3,
      strokewidth=2, label="Tile 2 footprint")
scatter!(axis, fused_transcripts; color=:black, markersize=13,
         label="fused transcript")
text!(axis, 6.6, 4.1; text="origin ambiguous", color=:purple,
      align=(:center, :bottom), fontsize=13)
axislegend(axis; position=:rt, framevisible=false)
xlims!(axis, -0.5, 13.5); ylims!(axis, -0.5, 9.5)
save("starmap-overlap.svg", figure)
Markdown.parse("![Overlapping registered STARmap tiles](starmap-overlap.svg)")
```

## Keep registered images tiled

A full bounding canvas wastes memory in gaps and silently needs a policy for overlapping
pixels. Keep each tile positioned instead:

```julia
import SpatialOmics as SO

function load_dapi_tiles(dataset, root)
    names = sources(dataset)
    tiles = SpatialImage[]
    for name in names
        attributes = source_attributes(dataset, name)
        raw = load(joinpath(root, attributes["dapi_path"]))
        projection = ndims(raw) == 3 ?
            dropdims(maximum(Float32.(raw); dims=3); dims=3) : Float32.(raw)
        x, y, _ = attributes["offset_xyz_px"]
        push!(tiles, SpatialImage(
            projection;
            coord_system="global_px",
            pixel_to_cs=SO.translation(x, y, "pixel", "global_px"),
        ))
    end
    SpatialRasterTiles(tiles, names)
end
```

`image!(axis, scaleminmax(tiles))` plots all pieces in the shared coordinate system.
The current native writer does not persist a `SpatialRasterTiles` collection. Until an
explicit tiled-raster export is implemented, keep the TIFF paths and offsets as the
recoverable representation rather than silently writing a dense, last-write-wins mosaic.

## Version the local reader

Treat the code above as a reader for a specific observed layout. In production, record a
reader version in dataset metadata and add small fixture tests for:

- required filenames and columns;
- coordinate conversion at tile corners;
- preservation of canvas coordinates and z values;
- source attributes and footprint IDs;
- overlapping tile behavior; and
- reconstruction of an equivalent vendor transcript table.

When the STARmap export changes, add a new method or explicit layout option. Do not
silently reinterpret an old cache; rebuild it with the appropriate reader version.
