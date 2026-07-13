
# ── Backing store ─────────────────────────────────────────────────────────────

"""
    BackingStore(; path=nothing)

Disk location for a dataset's Zarr storage, with ownership tracking.

When `path` is `nothing`, a temporary directory is created and owned by this
store (deleted automatically when the parent `SpatialDataset` is garbage
collected or `close`d). When `path` is supplied, the directory is used as-is
and the store is *not* owned — no automatic cleanup occurs.

# See also
[`SpatialDataset`](@ref), [`keep!`](@ref), [`with_dataset`](@ref)
"""
mutable struct BackingStore
    path    :: String
    owned   :: Bool
    handles :: Dict{String, Any}    # element name → open zarr group handle
end

function BackingStore(; path=nothing)
    if path === nothing
        p = mktempdir(; prefix="spatialomics_")
        owned = true
    else
        p = abspath(path)
        mkpath(p)
        owned = false
    end
    _init_zarr_root(p)
    BackingStore(p, owned, Dict{String,Any}())
end

function _init_zarr_root(path::String)
    zarr_json = joinpath(path, "zarr.json")
    isfile(zarr_json) && return
    open(zarr_json, "w") do io
        write(io, """{"zarr_format":3,"node_type":"group","attributes":{"spatialdata_attrs":{"version":"0.2.0"}}}""")
    end
    # Write spatialomics_meta.json so new stores are not mistaken for Python SpatialData format
    open(joinpath(path, "spatialomics_meta.json"), "w") do io
        write(io, """{"coord_systems":[]}""")
    end
end

function _cleanup!(bs::BackingStore)
    bs.owned || return
    close.(values(bs.handles))
    rm(bs.path; recursive=true, force=true)
end

# ── Dataset ───────────────────────────────────────────────────────────────────

"""
    SpatialDataset(; path=nothing, metadata=Dict())

Root container for a spatial omics experiment.

Holds named collections of spatial elements (`SpatialPoints`, `SpatialShapes`,
`SpatialImage`, `SpatialLabels`), a graph of `CoordinateSystem` nodes connected
by `AbstractTransformation` edges, named `SpatialRelation` objects, and free-form
metadata. Follows the [SpatialData specification](https://spatialdata.scverse.org/).

All data is backed by a `BackingStore` Zarr directory. When `path` is `nothing`,
a temporary directory is used and cleaned up automatically. Supply `path` to
write directly to a persistent location. Every `setindex!` call writes the
element to disk immediately — the dataset is always on disk.

```julia
ds = SpatialDataset()                          # temp-backed
ds = SpatialDataset(path="/data/exp.zarr")     # persistent-backed
```

# See also
[`BackingStore`](@ref), [`with_dataset`](@ref), [`keep!`](@ref),
[`elements`](@ref), [`coord_systems`](@ref), [`relations`](@ref)
"""
# ── Backed metadata dict ──────────────────────────────────────────────────────

"""
    BackedMetadata

Dict-like container for dataset metadata that writes each entry to the backing
store immediately on assignment, keeping disk and memory in sync.  Accessed
as `ds.metadata`.

String and integer keys are converted to `String` automatically.
"""
mutable struct BackedMetadata <: AbstractDict{String, Any}
    data    :: Dict{String, Any}
    backing :: BackingStore
end

function Base.setindex!(bm::BackedMetadata, val, key::String)
    bm.data[key] = val
    _write_metadata_entry(bm.backing.path, key, val)   # defined in zarr_io.jl
    bm
end

function Base.delete!(bm::BackedMetadata, key::String)
    delete!(bm.data, key)
    p = joinpath(bm.backing.path, "metadata", key)
    isdir(p) && rm(p; recursive=true, force=true)
    bm
end

Base.getindex(bm::BackedMetadata, key::String) = bm.data[key]
Base.iterate(bm::BackedMetadata)               = iterate(bm.data)
Base.iterate(bm::BackedMetadata, state)        = iterate(bm.data, state)
Base.length(bm::BackedMetadata)                = length(bm.data)

# ── Dataset ───────────────────────────────────────────────────────────────────

mutable struct SpatialDataset
    elements      :: OrderedDict{String, Any}
    coord_systems :: OrderedDict{String, CoordinateSystem}
    transforms    :: Vector{AbstractTransformation}
    backing       :: BackingStore
    relations     :: Dict{String, Any}     # name → SpatialRelation
    metadata      :: BackedMetadata
end

function SpatialDataset(; path=nothing, metadata=Dict{String,Any}())
    bs = BackingStore(; path)
    ds = SpatialDataset(
        OrderedDict{String,Any}(),
        OrderedDict{String,CoordinateSystem}(),
        AbstractTransformation[],
        bs,
        Dict{String,Any}(),
        BackedMetadata(Dict{String,Any}(metadata), bs),
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
    _write_spatialomics_meta(ds, ds.backing.path)   # defined in zarr_io.jl
    ds
end

function Base.push!(ds::SpatialDataset, t::AbstractTransformation)
    push!(ds.transforms, t)
    _write_spatialomics_meta(ds, ds.backing.path)   # defined in zarr_io.jl
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
    _write_zarr_relation(ds.backing.path, name, rel)
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
    _write_zarr(ds.backing.path, name, el)
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
