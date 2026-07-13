# ── Format token ──────────────────────────────────────────────────────────────

"""
    SpatialDataZarr()

Format token for the native SpatialData OME-Zarr on-disk format.

Pass to `read` or `write!` to select this backend:

```julia
ds = read(SpatialDataZarr(), "/path/to/experiment.zarr")
write!(ds, "/path/to/output.zarr", SpatialDataZarr())
```

`read` auto-detects whether the Zarr store was written by Python's SpatialData
library or by this package and dispatches accordingly.

# See also
[`CosMx`](@ref), [`write!`](@ref)
"""
struct SpatialDataZarr end

# ── Low-level zarr helpers ─────────────────────────────────────────────────────

function _write_group_meta(path::String, attrs::AbstractDict=Dict{String,Any}())
    meta = Dict("zarr_format" => 3, "node_type" => "group", "attributes" => attrs)
    open(joinpath(path, "zarr.json"), "w") do io
        JSON.print(io, meta)
    end
end

function _write_zarr_array(grp_path::String, name::String, arr::AbstractVector{T}) where T
    rm(joinpath(grp_path, name); recursive=true, force=true)
    z = zcreate(T, length(arr); path=grp_path, name, zarr_format=3,
                chunks=(length(arr),), fill_value=zero(T))
    z[:] = arr
end

const _ZARR_CHUNK_SIZE = 512

function _write_zarr_array(grp_path::String, name::String, arr::AbstractMatrix{T}) where T
    rm(joinpath(grp_path, name); recursive=true, force=true)
    m, n = size(arr)
    z = zcreate(T, m, n; path=grp_path, name, zarr_format=3,
                chunks=(min(_ZARR_CHUNK_SIZE, m), min(_ZARR_CHUNK_SIZE, n)),
                fill_value=zero(T))
    z[:, :] = arr
end

function _write_zarr_array(grp_path::String, name::String, arr::AbstractArray{T, 3}) where T
    rm(joinpath(grp_path, name); recursive=true, force=true)
    a, b, c = size(arr)
    z = zcreate(T, a, b, c; path=grp_path, name, zarr_format=3,
                chunks=(a, min(_ZARR_CHUNK_SIZE, b), min(_ZARR_CHUNK_SIZE, c)),
                fill_value=zero(T))
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
        JSON.print(io, pts.feature_codebook)
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
        convert(Vector{String}, JSON.parse(read(cb_path, String))) : String[]

    n = size(coords_mat, 1)
    coords = [Point2f(coords_mat[i, 1], coords_mat[i, 2]) for i in 1:n]

    meta = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    cs   = meta["attributes"]["_spatialdata_attrs"]["coord_system"]

    SpatialPoints{Float32}(coords, feature_id, codebook, instance_id, nothing, cs, nothing)
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

    meta = JSON.parse(read(joinpath(grp, "zarr.json"), String))
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
        mat  = SMatrix{3,3,Float64}(Float64(rows[i][j]) for i in 1:3, j in 1:3)
        Affine(mat, d["src"], d["dst"])
    else
        Identity(d["src"], d["dst"])
    end
end

# ── Write SpatialImage ─────────────────────────────────────────────────────────

function _write_zarr(root::String, name::String, img::SpatialImage{T}) where T
    grp = joinpath(root, "images", name)
    isdir(joinpath(grp, "data")) && return   # already written at this path — skip
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
    meta    = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    attrs   = meta["attributes"]["_spatialdata_attrs"]
    ax      = Tuple(Symbol.(attrs["axes"]))
    names   = String.(get(attrs, "channel_names", String[]))
    cs      = attrs["coord_system"]
    p2cs    = _transform_from_dict(attrs["pixel_to_cs"])
    N       = length(ax)

    img = SpatialImage(zopen(joinpath(grp, "data"), "r"; zarr_format=3);
                       axes=ax, channel_names=names, coord_system=cs, pixel_to_cs=p2cs)
    i = 1
    while isdir(joinpath(grp, "level$i"))
        push!(img.pyramid, zopen(joinpath(grp, "level$i"), "r"; zarr_format=3))
        i += 1
    end
    img
end

# ── Write SpatialLabels ────────────────────────────────────────────────────────

function _write_zarr(root::String, name::String, lbl::SpatialLabels{T}) where T
    grp = joinpath(root, "labels", name)
    mkpath(grp)
    _write_group_meta(grp, Dict(
        "_spatialdata_attrs" => Dict(
            "type"         => "labels",
            "axes"         => collect(string.(lbl.axes)),
            "coord_system" => lbl.coord_system,
            "pixel_to_cs"  => _transform_to_dict(lbl.pixel_to_cs))))
    _write_zarr_array(grp, "data", lbl.data)
    open(joinpath(grp, "instance_map.json"), "w") do io
        JSON.print(io, Dict(string(k) => v for (k, v) in lbl.instance_map))
    end
end

# ── Read SpatialLabels ─────────────────────────────────────────────────────────

function _read_labels_zarr(grp::String)
    meta  = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    attrs = meta["attributes"]["_spatialdata_attrs"]
    ax    = Tuple(Symbol.(attrs["axes"]))
    cs    = attrs["coord_system"]
    p2cs  = _transform_from_dict(attrs["pixel_to_cs"])
    N     = length(ax)

    data = zopen(joinpath(grp, "data"), "r"; zarr_format=3)
    T    = eltype(data)

    imap_path = joinpath(grp, "instance_map.json")
    instance_map = if isfile(imap_path)
        d = JSON.parse(read(imap_path, String))
        Dict{T, Int32}(parse(T, string(k)) => Int32(v) for (k, v) in d)
    else
        Dict{T, Int32}()
    end

    SpatialLabels(data; axes=ax, instance_map, coord_system=cs, pixel_to_cs=p2cs)
end

# ── JSON metadata helpers ──────────────────────────────────────────────────────

function _read_json_table(path::String)
    isfile(path) || return NamedTuple()
    d = JSON.parse(read(path, String))
    isempty(d) && return NamedTuple()
    cols = Tuple(Symbol.(keys(d)))
    vals = Tuple(collect(d[k]) for k in keys(d))
    NamedTuple{cols}(vals)
end

# ── Write SpatialRelation ──────────────────────────────────────────────────────

function _kind_meta(kind::Membership{strict}) where strict
    Dict("kind" => "Membership", "strict" => strict)
end
_kind_meta(::Expression) = Dict("kind" => "Expression")

_write_zarr(::String, ::String, ::Any) = nothing  # SpatialTable and future types not yet serialized

function _write_zarr_relation(root::String, name::String, rel::SpatialRelation)
    grp = joinpath(root, "relations", name)
    mkpath(grp)
    meta = merge(_kind_meta(rel.kind),
                 Dict("src" => rel.src,
                      "dst" => rel.dst === nothing ? nothing : rel.dst))
    _write_group_meta(grp, Dict("_spatialdata_attrs" => meta))
    _write_zarr_array(grp, "src_ids", rel.src_ids)
    isempty(rel.dst_ids) || _write_zarr_array(grp, "dst_ids", rel.dst_ids)
    rel.weights !== nothing && _write_zarr_array(grp, "weights", Matrix{Float32}(rel.weights))
    open(joinpath(grp, "obs.json"), "w") do io
        JSON.print(io, Dict(string(k) => collect(v) for (k, v) in pairs(Tables.columntable(rel.obs))))
    end
    open(joinpath(grp, "var.json"), "w") do io
        JSON.print(io, Dict(string(k) => collect(v) for (k, v) in pairs(Tables.columntable(rel.var))))
    end
end

# ── Read SpatialRelation ───────────────────────────────────────────────────────

function _kind_from_meta(d::AbstractDict)
    k = String(d["kind"])
    k == "Expression" && return Expression()
    strict = Bool(get(d, "strict", false))
    strict ? Membership{true}() : Membership{false}()
end

function _read_relation_zarr(grp::String)
    meta  = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    attrs = get(meta["attributes"], "_spatialdata_attrs", Dict{String,Any}())
    kind  = _kind_from_meta(attrs)
    src   = String(get(attrs, "src", ""))
    dst_v = get(attrs, "dst", nothing)
    dst   = dst_v === nothing ? nothing : String(dst_v)

    src_ids = Vector{Int32}(zopen(joinpath(grp, "src_ids"), "r"; zarr_format=3)[:])
    dst_ids = isfile(joinpath(grp, "dst_ids", "zarr.json")) ?
              Vector{Int32}(zopen(joinpath(grp, "dst_ids"), "r"; zarr_format=3)[:]) :
              Int32[]
    weights = isfile(joinpath(grp, "weights", "zarr.json")) ?
              zopen(joinpath(grp, "weights"), "r"; zarr_format=3)[:, :] :
              nothing
    obs = _read_json_table(joinpath(grp, "obs.json"))
    var = _read_json_table(joinpath(grp, "var.json"))
    SpatialRelation(src, dst, src_ids, dst_ids, weights, obs, var, kind)
end

# ── Dataset write ──────────────────────────────────────────────────────────────

function _write_metadata_entry(root::String, key::String, val::NamedTuple)
    grp = joinpath(root, "metadata", key)
    mkpath(grp)
    for (field, vec) in pairs(val)
        fname = string(field)
        if vec isa AbstractVector{<:Real}
            _write_zarr_array(grp, fname, collect(vec))
        elseif vec isa AbstractVector{<:AbstractString}
            codebook = unique(vec)
            idx_map  = Dict(s => Int32(i-1) for (i,s) in enumerate(codebook))
            ids      = Int32[idx_map[s] for s in vec]
            _write_zarr_array(grp, fname, ids)
            open(joinpath(grp, fname * "_codebook.json"), "w") do io
                JSON.print(io, codebook)
            end
        end
    end
    open(joinpath(grp, "type.json"), "w") do io
        JSON.print(io, Dict("type" => "named_tuple"))
    end
end

function _write_metadata_entry(root::String, key::String, val)
    grp = joinpath(root, "metadata", key)
    mkpath(grp)
    open(joinpath(grp, "value.json"), "w") do io
        JSON.print(io, val)
    end
end

function _write_all_metadata(ds::SpatialDataset, path::String)
    isempty(ds.metadata) && return
    for (key, val) in ds.metadata
        _write_metadata_entry(path, key, val)
    end
end

function _read_user_metadata!(ds::SpatialDataset, path::String)
    meta_dir = joinpath(path, "metadata")
    isdir(meta_dir) || return
    for key in readdir(meta_dir)
        grp = joinpath(meta_dir, key)
        isdir(grp) || continue
        type_path = joinpath(grp, "type.json")
        if isfile(type_path) && get(JSON.parse(read(type_path, String)), "type", "") == "named_tuple"
            fields = Symbol[]
            vecs   = AbstractVector[]
            for entry in readdir(grp)
                entry == "type.json"                && continue
                endswith(entry, "_codebook.json")   && continue
                !isdir(joinpath(grp, entry))        && continue
                codebook_path = joinpath(grp, entry * "_codebook.json")
                ids = zopen(joinpath(grp, entry), "r"; zarr_format=3)[:]
                if isfile(codebook_path)
                    codebook = convert(Vector{String},
                                       JSON.parse(read(codebook_path, String)))
                    push!(fields, Symbol(entry))
                    push!(vecs, [codebook[id+1] for id in ids])
                else
                    push!(fields, Symbol(entry))
                    push!(vecs, ids)
                end
            end
            isempty(fields) && continue
            ds.metadata.data[key] = NamedTuple{Tuple(fields)}(vecs)
        elseif isfile(joinpath(grp, "value.json"))
            ds.metadata.data[key] = JSON.parse(read(joinpath(grp, "value.json"), String))
        end
    end
end

function _write_spatialomics_meta(ds::SpatialDataset, path::String)
    open(joinpath(path, "spatialomics_meta.json"), "w") do io
        JSON.print(io, Dict(
            "coord_systems" => [
                Dict("name" => cs.name,
                     "axes"  => collect(string.(cs.axes)),
                     "units" => collect(cs.units))
                for cs in values(ds.coord_systems)],
            "transforms" => [_transform_to_dict(t) for t in ds.transforms]))
    end
end

function _write_dataset_zarr(ds::SpatialDataset, path::String)
    mkpath(path)
    _init_zarr_root(path)
    for subdir in ("points", "shapes", "images", "labels", "relations")
        mkpath(joinpath(path, subdir))
        _write_group_meta(joinpath(path, subdir))
    end
    for (name, el) in ds.elements
        _write_zarr(path, name, el)
    end
    for (name, rel) in ds.relations
        rel isa SpatialRelation && _write_zarr_relation(path, name, rel)
    end
    _write_spatialomics_meta(ds, path)
    _write_all_metadata(ds, path)
    path
end

function Base.write(ds::SpatialDataset, path::String, ::SpatialDataZarr)
    ds.backing.owned && @warn "Backing store is still at temp path \"$(ds.backing.path)\". " *
        "Call write!(ds, path, SpatialDataZarr()) to also update the dataset location."
    _write_dataset_zarr(ds, path)
end

"""
    write!(ds, path, SpatialDataZarr()) → ds

Write `ds` to the SpatialData OME-Zarr format at `path` and update the
dataset's backing store to point at the new location.

Unlike the non-mutating `write`, `write!` marks the backing store as permanent
(non-owned) so the directory is not deleted when `ds` is garbage collected.
Use this as the canonical "save" operation.

# See also
[`SpatialDataZarr`](@ref), [`keep!`](@ref)
"""
function write!(ds::SpatialDataset, path::String, ::SpatialDataZarr)
    _write_dataset_zarr(ds, path)
    ds.backing.path  = abspath(path)
    ds.backing.owned = false
    ds
end

# ── Dataset read ───────────────────────────────────────────────────────────────

function Base.read(::SpatialDataZarr, path::String) :: SpatialDataset
    isdir(path) || error("Path not found: $path")
    _is_python_spatialdata(path) && return _read_python_spatialdata(path)
    ds = SpatialDataset(; path)

    meta_path = joinpath(path, "spatialomics_meta.json")
    if isfile(meta_path)
        meta = JSON.parse(read(meta_path, String))
        for cs in get(meta, "coord_systems", [])
            ds.coord_systems[cs["name"]] = CoordinateSystem(cs["name"];
                                               axes  = Tuple(Symbol.(cs["axes"])),
                                               units = Tuple(String.(cs["units"])))
        end
        for t in get(meta, "transforms", [])
            push!(ds.transforms, _transform_from_dict(t))
        end
    end

    for (subdir, reader) in (("points", _read_points_zarr),
                              ("shapes", _read_shapes_zarr),
                              ("images", _read_image_zarr),
                              ("labels", _read_labels_zarr))
        for name in _zarr_element_names(path, subdir)
            ds.elements[name] = reader(joinpath(path, subdir, name))
        end
    end
    for name in _zarr_element_names(path, "relations")
        ds.relations[name] = _read_relation_zarr(joinpath(path, "relations", name))
    end
    _read_user_metadata!(ds, path)
    ds
end

function _zarr_element_names(root::String, kind::String)
    d = joinpath(root, kind)
    isdir(d) || return String[]
    filter(n -> isdir(joinpath(d, n)) && isfile(joinpath(d, n, "zarr.json")),
           readdir(d))
end

# ── Python SpatialData format — detection ─────────────────────────────────────

function _is_python_spatialdata(path::String)
    !isfile(joinpath(path, "spatialomics_meta.json")) &&
     isfile(joinpath(path, "zarr.json"))
end

# ── WKB polygon decoder (little-endian only) ───────────────────────────────────

function _wkb_polygon(bytes::AbstractVector{UInt8}) :: Union{Polygon, Nothing}
    io = IOBuffer(bytes)
    read(io, UInt8) == 0x01 || error("big-endian WKB not supported")
    geom_type = Int(read(io, UInt32))
    geom_type == 3 || return nothing   # only Polygon (type 3) supported; skip others
    n_rings = Int(read(io, UInt32))
    rings = map(1:n_rings) do _
        n_pts = Int(read(io, UInt32))
        [Point2f(read(io, Float64), read(io, Float64)) for _ in 1:n_pts]
    end
    length(rings) == 1 ? Polygon(rings[1]) : Polygon(rings[1], rings[2:end])
end

# ── OME-NGFF coord system helper ──────────────────────────────────────────────

function _ome_cs_name(attrs::AbstractDict)
    name = nothing
    for t in get(attrs, "coordinateTransformations", [])
        out = get(t, "output", nothing)
        out !== nothing && haskey(out, "name") && (name = String(out["name"]))
    end
    name !== nothing ? name : "global"
end

function _ome_xy_scale(attrs::AbstractDict)
    for t in get(attrs, "coordinateTransformations", [])
        get(t, "type", "") == "scale" || continue
        sc = get(t, "scale", nothing)
        sc === nothing && continue
        axes = [get(a, "name", "") for a in get(get(t, "input", Dict()), "axes", [])]
        if isempty(axes)
            length(sc) >= 2 && return Float64(sc[1]), Float64(sc[2])
        else
            xi = findfirst(==("x"), axes)
            yi = findfirst(==("y"), axes)
            xi !== nothing && yi !== nothing && return Float64(sc[xi]), Float64(sc[yi])
        end
    end
    1.0, 1.0
end

# ── Python OME-NGFF image reader ──────────────────────────────────────────────

function _read_ome_image_zarr_py(grp::String)
    meta = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    ome  = meta["attributes"]["ome"]
    ms   = ome["multiscales"][1]
    # Python zarr is C-order; Zarr.jl reverses dims → reverse axis labels too
    ax   = Tuple(reverse([Symbol(a["name"]) for a in ms["axes"]]))
    cs   = _ome_cs_name(ms)
    sx, sy = _ome_xy_scale(ms)
    p2cs = (sx ≈ 1.0 && sy ≈ 1.0) ? Identity("pixel", cs) : scaling(sx, sy, "pixel", cs)

    ch_names = String[]
    if haskey(get(ome, "omero", Dict()), "channels")
        ch_names = String[get(c, "label", "") for c in ome["omero"]["channels"]]
    end

    paths = String[d["path"] for d in ms["datasets"]]
    data  = zopen(joinpath(grp, paths[1]), "r"; zarr_format=3)
    img   = SpatialImage(data; axes=ax, channel_names=ch_names, coord_system=cs, pixel_to_cs=p2cs)
    for p in paths[2:end]
        push!(img.pyramid, zopen(joinpath(grp, p), "r"; zarr_format=3))
    end
    img
end

# ── Python OME-NGFF labels reader ─────────────────────────────────────────────

function _read_ome_labels_zarr_py(grp::String)
    meta = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    ome  = meta["attributes"]["ome"]
    ms   = ome["multiscales"][1]
    ax   = Tuple(reverse([Symbol(a["name"]) for a in ms["axes"]]))
    cs   = _ome_cs_name(ms)
    paths = String[d["path"] for d in ms["datasets"]]
    SpatialLabels(zopen(joinpath(grp, paths[1]), "r"; zarr_format=3); axes=ax, coord_system=cs)
end

# ── GeoParquet shapes reader ──────────────────────────────────────────────────

function _read_shapes_parquet(grp::String)
    pq_path = joinpath(grp, "shapes.parquet")
    meta    = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    cs      = _ome_cs_name(meta["attributes"])

    ds = Parquet2.Dataset(pq_path)
    rg = ds.row_groups[1]
    geoms_wkb = Parquet2.load(rg, 1)   # geometry column (always first)
    raw_ids   = Parquet2.load(rg, 2)   # index column (always second)

    maybe = [_wkb_polygon(g) for g in geoms_wkb]
    valid  = findall(!isnothing, maybe)
    isempty(valid) && error("No polygon geometries found in $(basename(grp))/shapes.parquet")

    polys  = Polygon[maybe[i] for i in valid]
    ids    = Int32.(1:length(valid))
    id_map = Dict{Int32, String}(Int32(j) => string(raw_ids[valid[j]]) for j in eachindex(valid))

    sx, sy = _ome_xy_scale(meta["attributes"])
    if !(sx ≈ 1.0 && sy ≈ 1.0)
        scale_pts(pts) = [Point2f(p[1]*sx, p[2]*sy) for p in pts]
        polys = [Polygon(scale_pts(p.exterior), [scale_pts(r) for r in p.interiors]) for p in polys]
    end

    SpatialShapes(polys; instance_id=ids, coord_system=cs), id_map
end

# ── Parquet points reader ─────────────────────────────────────────────────────

# cell_id columns may contain non-numeric strings (e.g. "UNASSIGNED") — map those to 0
_to_inst_id(v::Integer)        = Int32(v)
_to_inst_id(v::AbstractString) = something(tryparse(Int32, v), Int32(0))
_to_inst_id(v)                 = Int32(0)

function _read_points_parquet(grp::String)
    pq_dir   = joinpath(grp, "points.parquet")
    meta     = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    attrs    = meta["attributes"]
    sd_attrs = get(attrs, "spatialdata_attrs", Dict())
    feat_key = String(get(sd_attrs, "feature_key", "feature_name"))
    inst_key = String(get(sd_attrs, "instance_key", "cell_id"))
    cs       = _ome_cs_name(attrs)

    all_x    = Float32[]
    all_y    = Float32[]
    all_feat = String[]
    all_inst = Int32[]

    for part in sort(filter(f -> endswith(f, ".parquet"), readdir(pq_dir)))
        ds_part   = Parquet2.Dataset(joinpath(pq_dir, part))
        pq_meta   = JSON.parse(Parquet2.metadata(ds_part)["pandas"])
        col_names = String[c["name"] for c in pq_meta["columns"]]
        x_idx     = findfirst(==("x"),        col_names)
        y_idx     = findfirst(==("y"),        col_names)
        feat_idx  = findfirst(==(feat_key),   col_names)
        inst_idx  = findfirst(==(inst_key),   col_names)

        for rg in ds_part.row_groups
            xs = Parquet2.load(rg, x_idx)
            ys = Parquet2.load(rg, y_idx)
            ft = feat_idx !== nothing ? Parquet2.load(rg, feat_idx) : fill("", length(xs))
            id = inst_idx !== nothing ? Parquet2.load(rg, inst_idx) : zeros(Int64, length(xs))
            append!(all_x,    Float32.(xs))
            append!(all_y,    Float32.(ys))
            append!(all_feat, String.(ft))
            append!(all_inst, _to_inst_id.(id))
        end
    end

    sx, sy = _ome_xy_scale(attrs)
    if !(sx ≈ 1.0 && sy ≈ 1.0)
        all_x .*= Float32(sx)
        all_y .*= Float32(sy)
    end

    codebook   = sort(unique(all_feat))
    feat_to_id = Dict(g => Int32(i) for (i, g) in enumerate(codebook))
    feat_ids   = Int32[feat_to_id[f] for f in all_feat]

    SpatialPoints([Point2f(all_x[i], all_y[i]) for i in eachindex(all_x)];
                  feature_id=feat_ids, feature_codebook=codebook,
                  instance_id=all_inst, coord_system=cs)
end

# ── zarr v3 vlen-utf8 string array reader ─────────────────────────────────────
# Zarr.jl v0.10 cannot read data_type "string" arrays.
# This handles the AnnData vlen-utf8 + optional zstd codec chain (single-chunk).
# Binary layout after decompression: uint32 count, then per string: uint32 length + UTF-8 bytes.

function _decode_vlen_utf8(raw::Vector{UInt8})
    n   = Int(only(reinterpret(UInt32, raw[1:4])))
    pos = 5
    strings = Vector{String}(undef, n)
    for i in 1:n
        len = Int(only(reinterpret(UInt32, raw[pos:pos+3])))
        strings[i] = String(copy(raw[pos+4 : pos+3+len]))
        pos += 4 + len
    end
    strings
end

function _read_zarr_string_array(path::String)
    meta   = JSON.parse(read(joinpath(path, "zarr.json"), String))
    codecs = [c["name"] for c in get(meta, "codecs", [])]
    "vlen-utf8" in codecs || error("expected vlen-utf8 codec at $path; got $codecs")
    chunk_path = joinpath(path, "c", "0")
    isfile(chunk_path) || error("expected single chunk at $chunk_path; multi-chunk string arrays not supported")
    raw = read(chunk_path)
    _decode_vlen_utf8("zstd" in codecs ? transcode(ZstdDecompressor, raw) : raw)
end

# ── AnnData zarr table reader ─────────────────────────────────────────────────

function _read_anndata_table_zarr(grp::String)
    meta         = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    attrs        = meta["attributes"]
    region_str   = String(get(attrs, "region", ""))
    region       = isempty(region_str) ? nothing : region_str
    region_key   = Symbol(get(attrs, "region_key",   "region"))
    instance_key = Symbol(get(attrs, "instance_key", "instance_id"))

    X_data   = Vector{Float32}(zopen(joinpath(grp, "X", "data"),   "r"; zarr_format=3)[:])
    X_indices= Vector{Int32}(  zopen(joinpath(grp, "X", "indices"), "r"; zarr_format=3)[:])
    X_indptr = Vector{Int32}(  zopen(joinpath(grp, "X", "indptr"),  "r"; zarr_format=3)[:])
    X_meta   = JSON.parse(read(joinpath(grp, "X", "zarr.json"), String))
    n_obs, n_var = Int.(X_meta["attributes"]["shape"])

    # Build sparse row-index vector from CSR indptr, then construct SparseMatrixCSC.
    # Avoids the n_obs × n_var dense allocation that causes OOM on large tables.
    nnz  = length(X_data)
    rows = Vector{Int32}(undef, nnz)
    for i in 1:n_obs
        for j in (Int(X_indptr[i])+1):Int(X_indptr[i+1])
            rows[j] = i
        end
    end
    X = sparse(rows, X_indices .+ Int32(1), X_data, n_obs, n_var)

    _zarr_dtype(path) = get(JSON.parse(read(joinpath(path, "zarr.json"), String)), "data_type", "")

    var_index_path = joinpath(grp, "var", "_index")
    vnames = if isdir(var_index_path)
        _zarr_dtype(var_index_path) == "string" ?
            _read_zarr_string_array(var_index_path) :
            String.(zopen(var_index_path, "r"; zarr_format=3)[:])
    else
        String[]
    end
    var_nt = isempty(vnames) ? NamedTuple() : NamedTuple{(:name,)}((vnames,))

    # obs instance IDs: integer arrays → src_ids directly.
    # String arrays (e.g. Xenium cell barcodes) → synthetic Int32 src_ids + obs.name.
    obs_ids_path = joinpath(grp, "obs", string(instance_key))
    src_ids, obs_nt = if isdir(obs_ids_path)
        if _zarr_dtype(obs_ids_path) == "string"
            strs = _read_zarr_string_array(obs_ids_path)
            Int32.(1:n_obs), NamedTuple{(:name,)}((strs,))
        else
            Vector{Int32}(zopen(obs_ids_path, "r"; zarr_format=3)[:]), NamedTuple()
        end
    else
        Int32.(1:n_obs), NamedTuple()
    end

    src_str = isnothing(region) ? "" : region
    SpatialRelation(Expression(), src_str, src_ids, X; obs=obs_nt, var=var_nt)
end

# ── Python SpatialData dataset reader ─────────────────────────────────────────

function _read_python_spatialdata(path::String)
    ds = SpatialDataset(; path)
    str_id_maps = Dict{String, Dict{Int32, String}}()

    for (kind, reader) in (("images", _read_ome_image_zarr_py),
                            ("labels", _read_ome_labels_zarr_py))
        for name in _zarr_element_names(path, kind)
            try
                ds.elements[name] = reader(joinpath(path, kind, name))
            catch e
                @warn "Could not read $kind \"$name\": $e"
            end
        end
    end

    for name in _zarr_element_names(path, "tables")
        try
            ds.relations[name] = _read_anndata_table_zarr(joinpath(path, "tables", name))
        catch e
            @warn "Could not read table \"$name\": $e"
        end
    end

    for name in _zarr_element_names(path, "shapes")
        try
            el, id_map = _read_shapes_parquet(joinpath(path, "shapes", name))
            ds.elements[name] = el
            !isempty(id_map) && (str_id_maps[name] = id_map)
        catch e
            @warn "Could not read shapes \"$name\": $e"
        end
    end

    for name in _zarr_element_names(path, "points")
        try
            ds.elements[name] = _read_points_parquet(joinpath(path, "points", name))
        catch e
            @warn "Could not read points \"$name\": $e"
        end
    end

    !isempty(str_id_maps) && (ds.metadata["_instance_id_str_map"] = str_id_maps)
    # Register coord systems inferred from element metadata (Python format lacks explicit registry).
    # Bypass push!(ds, ...) to avoid writing spatialomics_meta.json into the Python SpatialData store.
    for (_, el) in ds.elements
        cs = coord_system(el)
        isempty(cs) && continue
        haskey(ds.coord_systems, cs) || (ds.coord_systems[cs] = CoordinateSystem(cs))
    end
    ds
end


# ── CosMx flatFiles reader ────────────────────────────────────────────────────

"""
    CosMx(; morphology_dir=nothing)

Format token for loading a CosMx SMI raw flat-file export.

Pass to `read` to load a CosMx export directory. If `morphology_dir` points to
a Morphology2D tile directory, tissue images are stitched and included as a
`SpatialImage`; otherwise only transcripts and cell boundaries are loaded.

```julia
ds = read(CosMx(), "/path/to/cosmx_export/")
ds = read(CosMx(morphology_dir="/path/to/Morphology2D"), "/path/to/cosmx_export/")
```

Each field-of-view (FOV) is registered as a separate `CoordinateSystem`;
use `coord_systems(ds)` and `transform(ds, fov_cs, "global")` to navigate
between spaces.

# See also
[`SpatialDataZarr`](@ref)
"""
Base.@kwdef struct CosMx
    morphology_dir :: Union{String, Nothing} = nothing
end

# Find the Morphology2D directory up to 4 levels below base_dir.
function _find_morphology2d(base_dir::String) :: Union{String, Nothing}
    isdir(base_dir) || return nothing
    function _search(d::String, depth::Int) :: Union{String, Nothing}
        depth == 0 && return nothing
        for entry in readdir(d; join=true)
            isdir(entry) || continue
            startswith(basename(entry), ".") && continue
            basename(entry) == "Morphology2D" && return entry
            found = _search(entry, depth - 1)
            found !== nothing && return found
        end
        nothing
    end
    _search(base_dir, 4)
end

# BFS down from base_dir (same root as _find_morphology2d) for the channel dict.
# Returns BiologicalTarget names in file order, or nothing if not found.
function _cosmx_channel_names(base_dir::String) :: Union{Vector{String}, Nothing}
    function _search(d::String, depth::Int) :: Union{Vector{String}, Nothing}
        depth == 0 && return nothing
        for entry in readdir(d; join=true)
            startswith(basename(entry), ".") && continue
            if isfile(entry) && basename(entry) == "Morphology_ChannelID_Dictionary.txt"
                lines = filter(!isempty, readlines(entry))
                length(lines) < 2 && return nothing
                return [split(l, '\t')[2] for l in lines[2:end]]
            end
            isdir(entry) && (r = _search(entry, depth - 1)) !== nothing && return r
        end
        nothing
    end
    isdir(base_dir) ? _search(base_dir, 6) : nothing
end

# Parse FOV ID from CosMx Morphology2D TIFF filenames: "*_F00001.TIF"
function _morphology_fov_tiffs(morph2d_dir::String) :: Dict{Int, String}
    tiff_map = Dict{Int, String}()
    for f in readdir(morph2d_dir)
        startswith(f, ".") && continue
        m = match(r"_F(\d+)\.(TIF|TIFF)$"i, f)
        m !== nothing && (tiff_map[parse(Int, m[1])] = joinpath(morph2d_dir, f))
    end
    tiff_map
end

# Function barrier: raw TIFF goes out of scope on return, making it GC-eligible
# before the next FOV is loaded. Avoids accumulating 40×168 MB in memory.
function _write_fov_morph!(z, tif_path, fov_x_val, fov_y_val, xmin, ymin, fov_w, fov_h, n_ch)
    raw = load(tif_path)
    cx = (round(Int, fov_x_val) - xmin + 1):(round(Int, fov_x_val) - xmin + fov_w)
    cy_bot = round(Int, fov_y_val) - fov_h + 1 - ymin + 1
    cy = cy_bot:(cy_bot + fov_h - 1)
    for ch in 1:n_ch
        page = ndims(raw) == 3 ? raw[:, :, ch] : raw
        z[cx, cy, ch] = Float32.(page[end:-1:1, :])'
    end
    nothing
end

# Stitch per-FOV Morphology2D TIFFs into a zarr at cache_path.
# Julia array layout: (canvas_w, canvas_h, n_ch) = (x, y, c).
function _stitch_morphology_to_zarr(morph2d_dir::String,
                                     fov_x::Dict{Int,Float64},
                                     fov_y::Dict{Int,Float64},
                                     cache_path::String;
                                     channel_names::Union{Vector{String}, Nothing}=nothing)
    tiff_map = _morphology_fov_tiffs(morph2d_dir)
    isempty(tiff_map) && error("No TIFF files found in $morph2d_dir")

    fov_ids = sort(collect(intersect(keys(tiff_map), keys(fov_x))))
    isempty(fov_ids) && error("No FOVs with both TIFF and position data")

    @info "Stitching $(length(fov_ids)) FOV TIFFs → $cache_path"

    # Load first TIFF to determine dimensions and channel count
    first_raw = load(tiff_map[fov_ids[1]])
    fov_h, fov_w = size(first_raw, 1), size(first_raw, 2)
    n_ch = ndims(first_raw) == 3 ? size(first_raw, 3) : 1

    # Canvas bounds in global_px (integer pixel coords)
    xmin = round(Int, minimum(fov_x[f] for f in fov_ids))
    xmax = round(Int, maximum(fov_x[f] for f in fov_ids)) + fov_w - 1
    # fov_y[f] = top edge of FOV (max global y); bottom = fov_y[f] - fov_h + 1
    ymin = round(Int, minimum(fov_y[f] for f in fov_ids)) - fov_h + 1
    ymax = round(Int, maximum(fov_y[f] for f in fov_ids))
    canvas_w = xmax - xmin + 1
    canvas_h = ymax - ymin + 1

    # Create zarr store and group metadata
    mkpath(cache_path)
    ch_names = if channel_names !== nothing && length(channel_names) == n_ch
        channel_names
    else
        channel_names !== nothing &&
            @warn "Channel name count ($(length(channel_names))) ≠ TIFF channels ($n_ch); using defaults"
        ["channel_$i" for i in 1:n_ch]
    end
    p2cs = translation(Float64(xmin), Float64(ymin), "pixel", "global_px")
    _write_group_meta(cache_path, Dict(
        "_spatialdata_attrs" => Dict(
            "type"          => "image",
            "axes"          => ["x", "y", "c"],
            "channel_names" => ch_names,
            "coord_system"  => "global_px",
            "pixel_to_cs"   => _transform_to_dict(p2cs))))

    # Julia array (canvas_w, canvas_h, n_ch) = (x, y, c); one channel per chunk
    # so each z[cx, cy, ch] write touches exactly one chunk with no read-modify-write.
    chunk_y = min(fov_h, canvas_h)
    chunk_x = min(fov_w, canvas_w)
    z0_path = joinpath(cache_path, "data")
    rm(z0_path; recursive=true, force=true)
    mkpath(z0_path)
    # Zarr.jl v3 reverses shape/chunks on read: zarr.json [A,B,C] → Julia size (C,B,A).
    # We want Julia (x,y,c) = (canvas_w,canvas_h,n_ch), so zarr.json = [n_ch,canvas_h,canvas_w].
    # Chunk reversal: Julia chunks (chunk_x,chunk_y,1) → zarr.json [1,chunk_y,chunk_x].
    open(joinpath(z0_path, "zarr.json"), "w") do io
        JSON.print(io, Dict(
            "zarr_format"        => 3,
            "node_type"          => "array",
            "shape"              => [n_ch, canvas_h, canvas_w],
            "data_type"          => "float32",
            "chunk_grid"         => Dict("name" => "regular",
                                         "configuration" => Dict("chunk_shape" => [1, chunk_y, chunk_x])),
            "chunk_key_encoding" => Dict("name" => "default",
                                         "configuration" => Dict("separator" => "/")),
            "fill_value"         => 0.0,
            "codecs"             => [Dict("name" => "bytes",
                                          "configuration" => Dict("endian" => "little"))],
            "attributes"         => Dict{String,Any}()))
    end
    z = zopen(z0_path, "w"; zarr_format=3)

    first_raw = nothing  # release before loop

    for f in fov_ids
        _write_fov_morph!(z, tiff_map[f], fov_x[f], fov_y[f], xmin, ymin, fov_w, fov_h, n_ch)
        GC.gc(false)  # free the loaded TIFF before the next iteration
    end
    @info "Morphology stitched to $cache_path"
    nothing
end

function _read_cosmx_morphology(zarr_path::String) :: SpatialImage{Float32}
    meta   = JSON.parse(read(joinpath(zarr_path, "zarr.json"), String))
    attrs  = meta["attributes"]["_spatialdata_attrs"]
    ax     = Tuple(Symbol.(attrs["axes"]))
    names  = String.(get(attrs, "channel_names", String[]))
    cs     = attrs["coord_system"]
    p2cs   = _transform_from_dict(attrs["pixel_to_cs"])

    data = zopen(joinpath(zarr_path, "data"), "r"; zarr_format=3)
    img  = SpatialImage(data; axes=ax, channel_names=names, coord_system=cs, pixel_to_cs=p2cs)
    i = 1
    while isdir(joinpath(zarr_path, "level$i"))
        push!(img.pyramid, zopen(joinpath(zarr_path, "level$i"), "r"; zarr_format=3))
        i += 1
    end
    img
end

# Returns the run directory containing tx_file, fov_positions, polygons CSVs.
# Accepts: the run dir directly, or a parent dir containing a single run dir.
function _cosmx_run_dir(path::String)
    isdir(path) || error("CosMx path not found: $path")
    # if path itself contains the expected files, use it directly
    if any(endswith(f, "_tx_file.csv.gz") for f in readdir(path))
        return path
    end
    # one level down
    for sub in readdir(path; join=true)
        isdir(sub) || continue
        if any(endswith(f, "_tx_file.csv.gz") for f in readdir(sub))
            return sub
        end
    end
    error("No CosMx tx_file.csv.gz found under $path")
end

function _cosmx_find(dir::String, suffix::String)
    hits = filter(f -> endswith(f, suffix) && !startswith(f, "._"), readdir(dir))
    isempty(hits) && error("No file with suffix $suffix in $dir")
    joinpath(dir, first(hits))
end

function _gz_csv(path::String)
    CSV.File(path)   # CSV.jl detects .gz extension and decompresses automatically
end

function Base.read(fmt::CosMx, path::String;
                   cache::Union{String, Nothing}=nothing) :: SpatialDataset
    run_dir = _cosmx_run_dir(path)

    # ── FOV positions ─────────────────────────────────────────────────────────
    fov_tbl = _gz_csv(_cosmx_find(run_dir, "_fov_positions_file.csv.gz"))
    fov_x   = Dict{Int, Float64}()
    fov_y   = Dict{Int, Float64}()
    for row in fov_tbl
        fov_x[Int(row.FOV)] = Float64(row.x_global_px)
        fov_y[Int(row.FOV)] = Float64(row.y_global_px)
    end
    fov_ids = sort(collect(keys(fov_x)))

    # ── Transcripts ───────────────────────────────────────────────────────────
    tx_tbl = _gz_csv(_cosmx_find(run_dir, "_tx_file.csv.gz"))

    all_x        = Float32[]
    all_y        = Float32[]
    all_feat     = String[]
    ann_fov      = Int32[]
    ann_z        = Float32[]
    ann_comp     = String[]
    ann_cell_id  = Int32[]   # raw per-FOV cell_ID; remapped to global after cells are read
    fov_max_local = Dict{Int, Tuple{Float32, Float32}}()

    for row in tx_tbl
        push!(all_x,       Float32(row.x_global_px))
        push!(all_y,       Float32(row.y_global_px))
        push!(all_feat,    String(row.target))
        push!(ann_fov,     Int32(row.fov))
        push!(ann_z,       row.z === missing ? Float32(0) : Float32(row.z))
        push!(ann_comp,    row.CellComp === missing ? "" : String(row.CellComp))
        push!(ann_cell_id, row.cell_ID === missing ? Int32(0) : Int32(row.cell_ID))
        f  = Int(row.fov)
        lx = Float32(row.x_global_px) - Float32(fov_x[f])
        ly = Float32(fov_y[f]) - Float32(row.y_global_px)
        prev = get(fov_max_local, f, (0f0, 0f0))
        fov_max_local[f] = (max(prev[1], lx), max(prev[2], ly))
    end

    fov_w = isempty(fov_max_local) ? 4256f0 : maximum(first, values(fov_max_local))
    fov_h = isempty(fov_max_local) ? 4256f0 : maximum(last,  values(fov_max_local))

    codebook   = sort(unique(all_feat))
    feat_to_id = Dict(g => Int32(i) for (i, g) in enumerate(codebook))
    feat_ids   = Int32[feat_to_id[f] for f in all_feat]

    # ── Cell polygons ─────────────────────────────────────────────────────────
    poly_tbl = _gz_csv(_cosmx_find(run_dir, "-polygons.csv.gz"))

    cell_keys  = Pair{Int,Int}[]
    key_set    = Set{Pair{Int,Int}}()
    vert_map   = Dict{Pair{Int,Int}, Vector{Point2f}}()

    for row in poly_tbl
        k = Int(row.fov) => Int(row.cellID)
        if k ∉ key_set
            push!(key_set, k)
            push!(cell_keys, k)
            vert_map[k] = Point2f[]
        end
        push!(vert_map[k], Point2f(Float32(row.x_global_px), Float32(row.y_global_px)))
    end

    # Sequential global IDs: (fov, cellID) → Int32 (1-based, unique across all FOVs)
    global_id = Dict{Pair{Int,Int}, Int32}(k => Int32(i) for (i, k) in enumerate(cell_keys))

    polys = Polygon[]
    inst  = Int32[]
    for (i, k) in enumerate(cell_keys)
        verts = vert_map[k]
        first(verts) ≈ last(verts) || push!(verts, verts[1])
        push!(polys, Polygon(verts))
        push!(inst, Int32(i))
    end

    cells = SpatialShapes(polys; instance_id=inst, coord_system="global_px")

    # Remap transcript instance_ids now that global_id map is available
    all_inst = Int32[ann_cell_id[i] == Int32(0) ? Int32(0) :
                     get(global_id, Int(ann_fov[i]) => Int(ann_cell_id[i]), Int32(0))
                     for i in eachindex(ann_cell_id)]

    transcripts = SpatialPoints(
        [Point2f(all_x[i], all_y[i]) for i in eachindex(all_x)];
        feature_id       = feat_ids,
        feature_codebook = codebook,
        instance_id      = all_inst,
        coord_system     = "global_px")

    # ── Assemble dataset ──────────────────────────────────────────────────────
    ds = cache !== nothing ? SpatialDataset(; path=abspath(cache)) : SpatialDataset()
    push!(ds, CoordinateSystem("global_px"; axes=(:x, :y), units=("px", "px")))

    for f in fov_ids
        cs_name = "fov_$(f)_px"
        push!(ds, CoordinateSystem(cs_name; axes=(:x, :y), units=("px", "px")))
        # local→global: x_g = x_l + fov_x[f], y_g = fov_y[f] - y_l
        # SMatrix fills column-by-column:
        #   col1=[1,0,0], col2=[0,-1,0], col3=[fov_x,fov_y,1]
        #   → rows: [1,0,fov_x], [0,-1,fov_y], [0,0,1]
        mat = SMatrix{3,3,Float64}(1, 0, 0,
                                   0,-1, 0,
                                   fov_x[f], fov_y[f], 1)
        push!(ds, Affine(mat, cs_name, "global_px"))
    end

    ds["transcripts"] = transcripts
    ds["cells"]       = cells
    ds["fovs"]        = SpatialShapes(
        [let ox = Float32(fov_x[f]), oy = Float32(fov_y[f])
             Polygon([Point2f(ox,         oy - fov_h),
                      Point2f(ox + fov_w, oy - fov_h),
                      Point2f(ox + fov_w, oy),
                      Point2f(ox,         oy),
                      Point2f(ox,         oy - fov_h)])
         end for f in fov_ids];
        instance_id = Int32.(fov_ids), coord_system = "global_px")

    ds.metadata["transcripts_annotations"] = (
        fov      = ann_fov,
        z        = ann_z,
        CellComp = ann_comp)

    if fmt.morphology_dir !== nothing
        morph2d = _find_morphology2d(fmt.morphology_dir)
        if morph2d === nothing
            @warn "No Morphology2D directory found under $(fmt.morphology_dir); skipping morphology"
        else
            morph_zarr = joinpath(ds.backing.path, "images", "morphology")
            isdir(morph_zarr) || _stitch_morphology_to_zarr(morph2d, fov_x, fov_y, morph_zarr;
                channel_names=_cosmx_channel_names(fmt.morphology_dir))
            ds["morphology"] = _read_cosmx_morphology(morph_zarr)
        end
    end

    cache !== nothing && _write_spatialomics_meta(ds, ds.backing.path)

    ds
end
