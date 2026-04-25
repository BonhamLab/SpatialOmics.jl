
# ── Backing store ─────────────────────────────────────────────────────────────

mutable struct BackingStore
    path            :: String
    owned           :: Bool
    spill_threshold :: Int                  # bytes; arrays larger than this spill to disk
    handles         :: Dict{String, Any}    # element name → open zarr group handle
end

function BackingStore(; path=nothing, spill_threshold=64_000_000)
    if path === nothing
        p = mktempdir(; prefix="spatialomics_")
        owned = true
    else
        p = abspath(path)
        mkpath(p)
        owned = false
    end
    _init_zarr_root(p)
    BackingStore(p, owned, spill_threshold, Dict{String,Any}())
end

function _init_zarr_root(path::String)
    # Write minimal SpatialData-compatible zarr root metadata
    zarr_json = joinpath(path, "zarr.json")
    isfile(zarr_json) && return
    open(zarr_json, "w") do io
        write(io, """{"zarr_format":3,"node_type":"group","attributes":{"spatialdata_attrs":{"version":"0.2.0"}}}""")
    end
end

function _cleanup!(bs::BackingStore)
    bs.owned || return
    try
        close.(values(bs.handles))
    catch
    end
    rm(bs.path; recursive=true, force=true)
end

# ── Dataset ───────────────────────────────────────────────────────────────────

mutable struct SpatialDataset
    elements      :: OrderedDict{String, Any}
    coord_systems :: OrderedDict{String, CoordinateSystem}
    transforms    :: Vector{AbstractTransformation}
    backing       :: BackingStore
    annotations   :: Dict{String, Any}     # element name → Tables.AbstractColumns
    metadata      :: Dict{String, Any}
end

function SpatialDataset(; path=nothing, spill_threshold=64_000_000, metadata=Dict{String,Any}())
    ds = SpatialDataset(
        OrderedDict{String,Any}(),
        OrderedDict{String,CoordinateSystem}(),
        AbstractTransformation[],
        BackingStore(; path, spill_threshold),
        Dict{String,Any}(),
        Dict{String,Any}(metadata),
    )
    finalizer(ds) do d
        d.backing.owned && _cleanup!(d.backing)
    end
    ds
end

# ── Lifecycle ─────────────────────────────────────────────────────────────────

function Base.close(ds::SpatialDataset)
    _cleanup!(ds.backing)
    ds.backing.owned = false   # prevent double-free in finalizer
    nothing
end

function keep!(ds::SpatialDataset, path::String=ds.backing.path)
    if path != ds.backing.path
        cp(ds.backing.path, path; force=true)
        ds.backing.owned && rm(ds.backing.path; recursive=true, force=true)
        ds.backing.path = abspath(path)
    end
    ds.backing.owned = false
    ds
end

function with_dataset(f::Function; path=nothing, kw...)
    ds = SpatialDataset(; path, kw...)
    try
        f(ds)
    finally
        close(ds)
    end
end

# ── Coordinate systems + transform API ────────────────────────────────────────

function Base.push!(ds::SpatialDataset, cs::CoordinateSystem)
    ds.coord_systems[cs.name] = cs
    ds
end

function Base.push!(ds::SpatialDataset, t::AbstractTransformation)
    push!(ds.transforms, t)
    ds
end

coord_systems(ds::SpatialDataset) = collect(keys(ds.coord_systems))

function transform(ds::SpatialDataset, src::String, dst::String)
    resolve_transform(ds.transforms, src, dst)
end

# ── Element attachment placeholder (implemented in elements.jl) ───────────────

function Base.setindex!(ds::SpatialDataset, el, name::String)
    _spill_element!(ds.backing, name, el)
    ds.elements[name] = el
    ds
end

Base.getindex(ds::SpatialDataset, name::String) = ds.elements[name]
Base.haskey(ds::SpatialDataset, name::String) = haskey(ds.elements, name)
Base.keys(ds::SpatialDataset) = keys(ds.elements)
