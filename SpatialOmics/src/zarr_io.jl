# ── Format token ──────────────────────────────────────────────────────────────

struct SpatialDataZarr end

# ── Low-level zarr helpers ─────────────────────────────────────────────────────

function _write_group_meta(path::String, attrs::AbstractDict=Dict{String,Any}())
    meta = Dict("zarr_format" => 3, "node_type" => "group", "attributes" => attrs)
    open(joinpath(path, "zarr.json"), "w") do io
        JSON3.write(io, meta)
    end
end

function _write_zarr_array(grp_path::String, name::String, arr::AbstractVector{T}) where T
    z = zcreate(T, length(arr); path=grp_path, name, zarr_format=3,
                chunks=(length(arr),), fill_value=zero(T))
    z[:] = arr
end

function _write_zarr_array(grp_path::String, name::String, arr::AbstractMatrix{T}) where T
    m, n = size(arr)
    z = zcreate(T, m, n; path=grp_path, name, zarr_format=3,
                chunks=(m, n), fill_value=zero(T))
    z[:, :] = arr
end

function _write_zarr_array(grp_path::String, name::String, arr::AbstractArray{T, 3}) where T
    a, b, c = size(arr)
    z = zcreate(T, a, b, c; path=grp_path, name, zarr_format=3,
                chunks=(a, b, c), fill_value=zero(T))
    z[:, :, :] = arr
end

function _read_zarr_array(grp_path::String, name::String)
    zopen(joinpath(grp_path, name), "r"; zarr_format=3)[fill(:, ndims(zopen(joinpath(grp_path, name), "r"; zarr_format=3)))...]
end

# ── Write SpatialPoints ────────────────────────────────────────────────────────

function _write_zarr(root::String, name::String, pts::SpatialPoints{T}) where T
    grp = joinpath(root, "points", name)
    mkpath(grp)
    _write_group_meta(grp, Dict(
        "_spatialdata_attrs" => Dict(
            "type" => "points",
            "axes" => [Dict("name" => "x", "type" => "space"),
                       Dict("name" => "y", "type" => "space")],
            "coord_system" => pts.coord_system)))

    n = length(pts)
    coords_mat = Matrix{Float32}(undef, n, 2)
    for i in 1:n
        coords_mat[i, 1] = pts.coords[i][1]
        coords_mat[i, 2] = pts.coords[i][2]
    end
    _write_zarr_array(grp, "coords",       coords_mat)
    _write_zarr_array(grp, "feature_id",   pts.feature_id)
    _write_zarr_array(grp, "instance_id",  pts.instance_id)

    open(joinpath(grp, "feature_codebook.json"), "w") do io
        JSON3.write(io, pts.feature_codebook)
    end
end

# ── Write SpatialShapes ────────────────────────────────────────────────────────

function _write_zarr(root::String, name::String, shp::SpatialShapes)
    grp = joinpath(root, "shapes", name)
    mkpath(grp)
    _write_group_meta(grp, Dict(
        "_spatialdata_attrs" => Dict(
            "type" => "shapes",
            "coord_system" => shp.coord_system)))

    _write_zarr_array(grp, "instance_id", shp.instance_id)

    # Ragged CSR layout: polygons → rings → points
    # poly_offsets[i] = 0-based index of first ring for polygon i (Julia 1-based)
    # ring_offsets[r] = 0-based index of first point for ring r (Julia 1-based)
    total_pts   = 0
    total_rings = 0
    for g in shp.geometries
        rings = GeoInterface.coordinates(g)
        total_rings += length(rings)
        for ring in rings
            total_pts += length(ring)
        end
    end

    geom_data     = Matrix{Float64}(undef, total_pts,   2)
    ring_offsets  = Vector{Int64}(undef,  total_rings + 1)
    poly_offsets  = Vector{Int64}(undef,  length(shp) + 1)

    pt_idx   = 0
    ring_idx = 0
    poly_offsets[1] = 0
    for (pi, g) in enumerate(shp.geometries)
        for ring in GeoInterface.coordinates(g)
            ring_offsets[ring_idx + 1] = pt_idx
            for pt in ring
                pt_idx += 1
                geom_data[pt_idx, 1] = Float64(pt[1])
                geom_data[pt_idx, 2] = Float64(pt[2])
            end
            ring_idx += 1
        end
        poly_offsets[pi + 1] = ring_idx
    end
    ring_offsets[end] = pt_idx

    _write_zarr_array(grp, "geom_data",    geom_data)
    _write_zarr_array(grp, "ring_offsets", ring_offsets)
    _write_zarr_array(grp, "poly_offsets", poly_offsets)
end

# ── Read SpatialPoints ─────────────────────────────────────────────────────────

function _read_points_zarr(grp::String) :: SpatialPoints{Float32}
    coords_mat   = zopen(joinpath(grp, "coords"),      "r"; zarr_format=3)[:, :]
    feature_id   = Vector{Int32}(zopen(joinpath(grp, "feature_id"),  "r"; zarr_format=3)[:])
    instance_id  = Vector{Int32}(zopen(joinpath(grp, "instance_id"), "r"; zarr_format=3)[:])

    cb_path = joinpath(grp, "feature_codebook.json")
    codebook = isfile(cb_path) ?
        JSON3.read(read(cb_path, String), Vector{String}) : String[]

    n = size(coords_mat, 1)
    coords = [Point2f(coords_mat[i, 1], coords_mat[i, 2]) for i in 1:n]

    meta = JSON3.read(read(joinpath(grp, "zarr.json"), String))
    cs   = meta["attributes"]["_spatialdata_attrs"]["coord_system"]

    SpatialPoints{Float32}(coords, feature_id, codebook, instance_id, cs)
end

# ── Read SpatialShapes ─────────────────────────────────────────────────────────

function _read_shapes_zarr(grp::String) :: SpatialShapes
    instance_id  = Vector{Int32}(zopen(joinpath(grp, "instance_id"),  "r"; zarr_format=3)[:])
    geom_data    = zopen(joinpath(grp, "geom_data"),    "r"; zarr_format=3)[:, :]
    ring_offsets = Vector{Int64}(zopen(joinpath(grp, "ring_offsets"), "r"; zarr_format=3)[:])
    poly_offsets = Vector{Int64}(zopen(joinpath(grp, "poly_offsets"), "r"; zarr_format=3)[:])

    n = length(instance_id)
    geometries = Vector{Polygon}(undef, n)
    for pi in 1:n
        r_start = poly_offsets[pi]     + 1    # 0-based offset → Julia 1-based start
        r_end   = poly_offsets[pi + 1]        # 0-based exclusive = Julia 1-based end
        rings   = Vector{Vector{Point2f}}(undef, r_end - r_start + 1)
        for (ri, r_idx) in enumerate(r_start:r_end)
            pt_start = ring_offsets[r_idx]     + 1
            pt_end   = ring_offsets[r_idx + 1]
            rings[ri] = [Point2f(geom_data[j, 1], geom_data[j, 2]) for j in pt_start:pt_end]
        end
        geometries[pi] = length(rings) == 1 ? Polygon(rings[1]) : Polygon(rings[1], rings[2:end])
    end

    meta = JSON3.read(read(joinpath(grp, "zarr.json"), String))
    cs   = meta["attributes"]["_spatialdata_attrs"]["coord_system"]

    SpatialShapes(geometries; instance_id, coord_system=cs)
end

# ── Transform serialization helpers ───────────────────────────────────────────

_transform_to_dict(t::Identity) =
    Dict("type" => "identity", "src" => t.src, "dst" => t.dst)

_transform_to_dict(t::Affine) =
    Dict("type" => "affine", "src" => t.src, "dst" => t.dst,
         "matrix" => [collect(t.matrix[i, :]) for i in 1:3])

_transform_to_dict(::AbstractTransformation) =
    Dict("type" => "identity", "src" => "", "dst" => "")

function _transform_from_dict(d)
    if d["type"] == "affine"
        rows = d["matrix"]
        mat  = SMatrix{3,3,Float64}(Float64(rows[j][i]) for i in 1:3, j in 1:3)
        Affine(mat, d["src"], d["dst"])
    else
        Identity(d["src"], d["dst"])
    end
end

# ── Write SpatialImage ─────────────────────────────────────────────────────────

function _write_zarr(root::String, name::String, img::SpatialImage{T}) where T
    grp = joinpath(root, "images", name)
    mkpath(grp)
    _write_group_meta(grp, Dict(
        "_spatialdata_attrs" => Dict(
            "type"          => "image",
            "axes"          => collect(string.(img.axes)),
            "channel_names" => img.channel_names,
            "coord_system"  => img.coord_system,
            "pixel_to_cs"   => _transform_to_dict(img.pixel_to_cs))))
    _write_zarr_array(grp, "data", img.data)
    for (i, level) in enumerate(img.pyramid)
        _write_zarr_array(grp, "level$(i)", level)
    end
end

# ── Read SpatialImage ──────────────────────────────────────────────────────────

function _read_image_zarr(grp::String)
    meta    = JSON3.read(read(joinpath(grp, "zarr.json"), String))
    attrs   = meta["attributes"]["_spatialdata_attrs"]
    ax      = Tuple(Symbol.(attrs["axes"]))
    names   = String.(get(attrs, "channel_names", String[]))
    cs      = attrs["coord_system"]
    p2cs    = _transform_from_dict(attrs["pixel_to_cs"])
    N       = length(ax)

    raw  = zopen(joinpath(grp, "data"), "r"; zarr_format=3)
    arr  = N == 2 ? raw[:, :] : raw[:, :, :]     # read into memory
    img  = SpatialImage(arr; axes=ax, channel_names=names, coord_system=cs, pixel_to_cs=p2cs)

    i = 1
    while isdir(joinpath(grp, "level$i"))
        raw_l = zopen(joinpath(grp, "level$i"), "r"; zarr_format=3)
        level = N == 2 ? raw_l[:, :] : raw_l[:, :, :]
        push!(img.pyramid, level)
        i += 1
    end
    img
end

# ── Dataset write ──────────────────────────────────────────────────────────────

function Base.write(ds::SpatialDataset, path::String, ::SpatialDataZarr)
    mkpath(path)
    _init_zarr_root(path)
    for subdir in ("points", "shapes", "images")
        mkpath(joinpath(path, subdir))
        _write_group_meta(joinpath(path, subdir))
    end
    for (name, el) in ds.elements
        if el isa SpatialPoints || el isa SpatialShapes || el isa SpatialImage
            _write_zarr(path, name, el)
        end
    end
    meta = Dict(
        "coord_systems" => [
            Dict("name" => cs.name,
                 "axes"  => collect(string.(cs.axes)),
                 "units" => collect(cs.units))
            for cs in values(ds.coord_systems)],
    )
    open(joinpath(path, "spatialomics_meta.json"), "w") do io
        JSON3.write(io, meta)
    end
    path
end

# ── Dataset read ───────────────────────────────────────────────────────────────

function Base.read(::SpatialDataZarr, path::String) :: SpatialDataset
    isdir(path) || error("Path not found: $path")
    ds = SpatialDataset(; path, spill_threshold=typemax(Int))

    meta_path = joinpath(path, "spatialomics_meta.json")
    if isfile(meta_path)
        meta = JSON3.read(read(meta_path, String))
        for cs in get(meta, "coord_systems", [])
            push!(ds, CoordinateSystem(cs["name"];
                                       axes  = Tuple(Symbol.(cs["axes"])),
                                       units = Tuple(String.(cs["units"]))))
        end
    end

    for (subdir, reader) in (("points", _read_points_zarr),
                              ("shapes", _read_shapes_zarr),
                              ("images", _read_image_zarr))
        for name in _zarr_element_names(path, subdir)
            ds.elements[name] = reader(joinpath(path, subdir, name))
        end
    end
    ds
end

function _zarr_element_names(root::String, kind::String)
    d = joinpath(root, kind)
    isdir(d) || return String[]
    filter(n -> isdir(joinpath(d, n)) && isfile(joinpath(d, n, "zarr.json")),
           readdir(d))
end

# ── Spill on attach ────────────────────────────────────────────────────────────

_element_bytes(pts::SpatialPoints{T}) where T =
    (2 * sizeof(T) + 2 * sizeof(Int32)) * length(pts)

_element_bytes(shp::SpatialShapes) =
    sizeof(Float64) * 4 * length(shp)   # bbox rows as rough lower bound

_element_bytes(::Any) = 0

function _spill_element!(bs::BackingStore, name::String, el)
    _element_bytes(el) < bs.spill_threshold && return
    (el isa SpatialPoints || el isa SpatialShapes) && _write_zarr(bs.path, name, el)
    nothing
end
