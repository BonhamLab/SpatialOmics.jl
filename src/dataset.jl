
# ── Backing store ─────────────────────────────────────────────────────────────

"""
    BackingStore(; path=nothing, spill_threshold=64_000_000)

Disk location for a dataset's Zarr storage, with ownership tracking.

When `path` is `nothing`, a temporary directory is created and owned by this
store (deleted automatically when the parent `SpatialDataset` is garbage
collected or `close`d). When `path` is supplied, the directory is used as-is
and the store is *not* owned — no automatic cleanup occurs.

`spill_threshold` (bytes) controls when large arrays are written to disk
immediately on element attachment rather than held in memory.

# See also
[`SpatialDataset`](@ref), [`keep!`](@ref), [`with_dataset`](@ref)
"""
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
    close.(values(bs.handles))
    rm(bs.path; recursive=true, force=true)
end

# ── Dataset ───────────────────────────────────────────────────────────────────

"""
    SpatialDataset(; path=nothing, spill_threshold=64_000_000, metadata=Dict())

Root container for a spatial omics experiment.

Holds named collections of spatial elements (`SpatialPoints`, `SpatialShapes`,
`SpatialImage`, `SpatialLabels`), a graph of `CoordinateSystem` nodes connected
by `AbstractTransformation` edges, named `SpatialRelation` objects, and free-form
metadata. Follows the [SpatialData specification](https://spatialdata.scverse.org/).

All data is backed by a `BackingStore` Zarr directory. When `path` is `nothing`,
a temporary directory is used and cleaned up automatically. Supply `path` to
write directly to a persistent location.

```julia
ds = SpatialDataset()                          # temp-backed
ds = SpatialDataset(path="/data/exp.zarr")     # persistent-backed
```

# See also
[`BackingStore`](@ref), [`with_dataset`](@ref), [`keep!`](@ref),
[`elements`](@ref), [`coord_systems`](@ref), [`relations`](@ref)
"""
mutable struct SpatialDataset
    elements      :: OrderedDict{String, Any}
    coord_systems :: OrderedDict{String, CoordinateSystem}
    transforms    :: Vector{AbstractTransformation}
    backing       :: BackingStore
    relations     :: Dict{String, Any}     # name → SpatialRelation
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

"""
    keep!(ds, path=ds.backing.path) → ds

Mark the dataset's backing store as permanent, preventing automatic cleanup.

If `path` differs from the current backing path, the store is copied there
first. After `keep!`, the dataset no longer owns its backing directory — it
will not be deleted when `ds` is garbage collected or `close`d.

# See also
[`with_dataset`](@ref), [`write!`](@ref)
"""
function keep!(ds::SpatialDataset, path::String=ds.backing.path)
    if path != ds.backing.path
        cp(ds.backing.path, path; force=true)
        ds.backing.owned && rm(ds.backing.path; recursive=true, force=true)
        ds.backing.path = abspath(path)
    end
    ds.backing.owned = false
    ds
end

"""
    with_dataset(f; path=nothing, kw...)

Open a dataset, run `f(ds)`, then close and clean up the backing store.

The dataset is always closed in a `finally` block, making this safe for
temporary analysis workflows that should not leave stale Zarr directories on disk.

```julia
result = with_dataset() do ds
    ds["cells"] = cells
    analyze(Expression(), transcripts, cells)
end
```

# See also
[`keep!`](@ref), [`SpatialDataset`](@ref)
"""
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

"""
    elements(ds) → OrderedDict{String, Any}

Return the ordered dictionary of all named spatial elements in `ds`.

Values are concrete element types (`SpatialPoints`, `SpatialShapes`,
`SpatialImage`, `SpatialLabels`). Use the typed accessors [`points`](@ref),
[`shapes`](@ref), [`images`](@ref), [`labels`](@ref) to retrieve a specific
element with type checking.
"""
elements(ds::SpatialDataset)      = ds.elements

"""
    coord_systems(ds) → Vector{String}

Return the names of all coordinate systems registered in `ds`.
"""
coord_systems(ds::SpatialDataset) = collect(keys(ds.coord_systems))

"""
    transform(ds, src, dst) → AbstractTransformation

Resolve a composed transformation from coordinate system `src` to `dst` using
the dataset's registered transforms. Delegates to [`resolve`](@ref).
"""
function transform(ds::SpatialDataset, src::String, dst::String)
    resolve(ds.transforms, src, dst)
end

# ── Element attachment placeholder (implemented in elements.jl) ───────────────

function Base.setindex!(ds::SpatialDataset, rel::SpatialRelation, name::String)
    _spill_relation!(ds.backing, name, rel)
    ds.relations[name] = rel
    ds
end

function Base.setindex!(ds::SpatialDataset, el, name::String)
    existing = _owning_dataset(el)
    if existing !== nothing && existing !== ds
        error("Element already attached to a different dataset. " *
              "Use `ds[\"$name\"] = copy(el)` to attach a detached copy. " *
              "Note: relations involving this element in the original dataset will not transfer.")
    end
    _spill_element!(ds.backing, name, el)
    _set_backref!(el, ds, name)
    ds.elements[name] = el
    ds
end

Base.getindex(ds::SpatialDataset, name::String) = ds.elements[name]
Base.haskey(ds::SpatialDataset, name::String) = haskey(ds.elements, name)
Base.keys(ds::SpatialDataset) = keys(ds.elements)

"""
    relations(ds) → Dict{String, SpatialRelation}
    relations(ds, name) → SpatialRelation

Return the dictionary of all named relations, or a specific relation by name.

# See also
[`SpatialRelation`](@ref), [`analyze`](@ref)
"""
function relations(ds::SpatialDataset)
    ds.relations
end

function relations(ds::SpatialDataset, name::String)
    haskey(ds.relations, name) ||
        error("No relation \"$name\". Available: $(join(keys(ds.relations), ", "))")
    ds.relations[name]
end
