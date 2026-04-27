# ── Format token ──────────────────────────────────────────────────────────────

struct SpatialDataZarr end

# ── Low-level zarr helpers ─────────────────────────────────────────────────────

function _write_group_meta(path::String, attrs::AbstractDict=Dict{String,Any}())
    meta = Dict("zarr_format" => 3, "node_type" => "group", "attributes" => attrs)
    open(joinpath(path, "zarr.json"), "w") do io
        JSON.print(io, meta)
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
    meta    = JSON.parse(read(joinpath(grp, "zarr.json"), String))
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

    raw  = zopen(joinpath(grp, "data"), "r"; zarr_format=3)
    data = N == 2 ? raw[:, :] : raw[:, :, :]
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

# ── Write SpatialTable ─────────────────────────────────────────────────────────

function _write_zarr(root::String, name::String, tbl::SpatialTable)
    grp = joinpath(root, "tables", name)
    mkpath(grp)
    _write_group_meta(grp, Dict(
        "_spatialdata_attrs" => Dict(
            "type"         => "table",
            "region"       => tbl.region === nothing ? "" : tbl.region,
            "region_key"   => string(tbl.region_key),
            "instance_key" => string(tbl.instance_key),
            "region_kind"  => string(tbl.region_kind))))
    _write_zarr_array(grp, "X", Matrix{Float32}(tbl.X))
    open(joinpath(grp, "obs.json"), "w") do io
        JSON.print(io, Dict(string(k) => collect(v) for (k, v) in pairs(tbl.obs)))
    end
    open(joinpath(grp, "var.json"), "w") do io
        JSON.print(io, Dict(string(k) => collect(v) for (k, v) in pairs(tbl.var)))
    end
end

# ── Read SpatialTable ──────────────────────────────────────────────────────────

function _read_json_table(path::String)
    isfile(path) || return NamedTuple()
    d = JSON.parse(read(path, String))
    isempty(d) && return NamedTuple()
    cols = Tuple(Symbol.(keys(d)))
    vals = Tuple(collect(d[k]) for k in keys(d))
    NamedTuple{cols}(vals)
end

function _read_table_zarr(grp::String)
    meta  = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    attrs = meta["attributes"]["_spatialdata_attrs"]
    region_str   = String(get(attrs, "region", ""))
    region       = isempty(region_str) ? nothing : region_str
    region_key   = Symbol(attrs["region_key"])
    instance_key = Symbol(attrs["instance_key"])
    region_kind  = Symbol(attrs["region_kind"])

    raw = zopen(joinpath(grp, "X"), "r"; zarr_format=3)
    X   = raw[:, :]
    obs = _read_json_table(joinpath(grp, "obs.json"))
    var = _read_json_table(joinpath(grp, "var.json"))

    SpatialTable(X; obs, var, region, region_key, instance_key, region_kind)
end

# ── Dataset write ──────────────────────────────────────────────────────────────

function Base.write(ds::SpatialDataset, path::String, ::SpatialDataZarr)
    mkpath(path)
    _init_zarr_root(path)
    for subdir in ("points", "shapes", "images", "labels", "tables")
        mkpath(joinpath(path, subdir))
        _write_group_meta(joinpath(path, subdir))
    end
    for (name, el) in ds.elements
        if el isa SpatialPoints || el isa SpatialShapes || el isa SpatialImage ||
           el isa SpatialLabels || el isa SpatialTable
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
        JSON.print(io, meta)
    end
    path
end

# ── Dataset read ───────────────────────────────────────────────────────────────

function Base.read(::SpatialDataZarr, path::String) :: SpatialDataset
    isdir(path) || error("Path not found: $path")
    _is_python_spatialdata(path) && return _read_python_spatialdata(path)
    ds = SpatialDataset(; path, spill_threshold=typemax(Int))

    meta_path = joinpath(path, "spatialomics_meta.json")
    if isfile(meta_path)
        meta = JSON.parse(read(meta_path, String))
        for cs in get(meta, "coord_systems", [])
            push!(ds, CoordinateSystem(cs["name"];
                                       axes  = Tuple(Symbol.(cs["axes"])),
                                       units = Tuple(String.(cs["units"]))))
        end
    end

    for (subdir, reader) in (("points", _read_points_zarr),
                              ("shapes", _read_shapes_zarr),
                              ("images", _read_image_zarr),
                              ("labels", _read_labels_zarr),
                              ("tables", _read_table_zarr))
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
    for t in get(attrs, "coordinateTransformations", [])
        out = get(t, "output", nothing)
        out !== nothing && haskey(out, "name") && return String(out["name"])
    end
    "global"
end

# ── Python OME-NGFF image reader ──────────────────────────────────────────────

function _read_ome_image_zarr_py(grp::String)
    meta = JSON.parse(read(joinpath(grp, "zarr.json"), String))
    ome  = meta["attributes"]["ome"]
    ms   = ome["multiscales"][1]
    # Python zarr is C-order; Zarr.jl reverses dims → reverse axis labels too
    ax   = Tuple(reverse([Symbol(a["name"]) for a in ms["axes"]]))
    cs   = _ome_cs_name(ms)

    ch_names = String[]
    if haskey(get(ome, "omero", Dict()), "channels")
        ch_names = String[get(c, "label", "") for c in ome["omero"]["channels"]]
    end

    paths = String[d["path"] for d in ms["datasets"]]
    data  = zopen(joinpath(grp, paths[1]), "r"; zarr_format=3)
    img   = SpatialImage(data; axes=ax, channel_names=ch_names, coord_system=cs)
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

    SpatialShapes(polys; instance_id=ids, coord_system=cs), id_map
end

# ── Parquet points reader ─────────────────────────────────────────────────────

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
            # cell_id may contain non-numeric strings (e.g. "UNASSIGNED") → 0
            append!(all_inst, [try Int32(v) catch; Int32(0) end for v in id])
        end
    end

    codebook   = sort(unique(all_feat))
    feat_to_id = Dict(g => Int32(i) for (i, g) in enumerate(codebook))
    feat_ids   = Int32[feat_to_id[f] for f in all_feat]

    SpatialPoints([Point2f(all_x[i], all_y[i]) for i in eachindex(all_x)];
                  feature_id=feat_ids, feature_codebook=codebook,
                  instance_id=all_inst, coord_system=cs)
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

    X = zeros(Float32, n_obs, n_var)
    for i in 1:n_obs
        for j in (Int(X_indptr[i])+1):Int(X_indptr[i+1])
            X[i, Int(X_indices[j])+1] += X_data[j]
        end
    end

    # Gene names from var/_index (zarr v3 string array)
    vnames = String[]
    try
        raw = zopen(joinpath(grp, "var", "_index"), "r"; zarr_format=3)
        vnames = String.(raw[:])
    catch
    end

    var_nt = isempty(vnames) ? NamedTuple() : (name = vnames,)
    SpatialTable(X; var=var_nt, region, region_key, instance_key, region_kind=:shapes)
end

# ── Python SpatialData dataset reader ─────────────────────────────────────────

function _read_python_spatialdata(path::String)
    ds = SpatialDataset(; path, spill_threshold=typemax(Int))
    str_id_maps = Dict{String, Dict{Int32, String}}()

    for (kind, reader) in (
            ("images",  _read_ome_image_zarr_py),
            ("labels",  _read_ome_labels_zarr_py),
            ("tables",  _read_anndata_table_zarr))
        for name in _zarr_element_names(path, kind)
            try
                ds.elements[name] = reader(joinpath(path, kind, name))
            catch e
                @warn "Could not read $kind \"$name\": $e"
            end
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
    ds
end

# ── Spill on attach ────────────────────────────────────────────────────────────

_element_bytes(pts::SpatialPoints{T}) where T =
    (2 * sizeof(T) + 2 * sizeof(Int32)) * length(pts)

_element_bytes(shp::SpatialShapes) =
    sizeof(Float64) * 4 * length(shp)   # bbox rows as rough lower bound

_element_bytes(tbl::SpatialTable) = sizeof(eltype(tbl.X)) * length(tbl.X)
_element_bytes(::Any) = 0

function _spill_element!(bs::BackingStore, name::String, el)
    _element_bytes(el) < bs.spill_threshold && return
    (el isa SpatialPoints || el isa SpatialShapes) && _write_zarr(bs.path, name, el)
    nothing
end

# ── CosMx flatFiles reader ────────────────────────────────────────────────────

struct CosMx end

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

function Base.read(::CosMx, path::String) :: SpatialDataset
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

    all_x    = Float32[]
    all_y    = Float32[]
    all_feat = String[]
    all_inst = Int32[]
    ann_fov  = Int32[]
    ann_z    = Float32[]
    ann_comp = String[]

    for row in tx_tbl
        push!(all_x,    Float32(row.x_global_px))
        push!(all_y,    Float32(row.y_global_px))
        push!(all_feat, String(row.target))
        push!(all_inst, row.cell_ID === missing ? Int32(0) : Int32(row.cell_ID))
        push!(ann_fov,  Int32(row.fov))
        push!(ann_z,    row.z === missing ? Float32(0) : Float32(row.z))
        push!(ann_comp, row.CellComp === missing ? "" : String(row.CellComp))
    end

    codebook   = sort(unique(all_feat))
    feat_to_id = Dict(g => Int32(i) for (i, g) in enumerate(codebook))
    feat_ids   = Int32[feat_to_id[f] for f in all_feat]

    transcripts = SpatialPoints(
        [Point2f(all_x[i], all_y[i]) for i in eachindex(all_x)];
        feature_id       = feat_ids,
        feature_codebook = codebook,
        instance_id      = all_inst,
        coord_system     = "global_px")

    # ── Cell polygons ─────────────────────────────────────────────────────────
    poly_tbl = _gz_csv(_cosmx_find(run_dir, "-polygons.csv.gz"))

    # Group vertices by (fov, cellID) preserving order
    cell_keys  = Pair{Int,Int}[]  # ordered unique (fov, cellID) pairs
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

    polys = Polygon[]
    inst  = Int32[]
    for (i, k) in enumerate(cell_keys)
        verts = vert_map[k]
        # ensure closed ring
        first(verts) ≈ last(verts) || push!(verts, verts[1])
        push!(polys, Polygon(verts))
        push!(inst, Int32(k.second))
    end

    cells = SpatialShapes(polys; instance_id=inst, coord_system="global_px")

    # ── Assemble dataset ──────────────────────────────────────────────────────
    ds = SpatialDataset()
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
        push!(ds.transforms, Affine(mat, cs_name, "global_px"))
    end

    ds["transcripts"] = transcripts
    ds["cells"]       = cells

    ds.metadata["transcripts_annotations"] = (
        fov      = ann_fov,
        z        = ann_z,
        CellComp = ann_comp)

    ds
end
