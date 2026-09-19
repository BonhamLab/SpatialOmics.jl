
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
    changes :: OrderedDict{Tuple{Symbol, String}, Symbol}
    closed  :: Bool
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
    BackingStore(
        p,
        owned,
        Dict{String,Any}(),
        OrderedDict{Tuple{Symbol,String},Symbol}(),
        false,
    )
end

function _init_zarr_root(path::String)
    zarr_json = joinpath(path, "zarr.json")
    isfile(zarr_json) && return
    open(zarr_json, "w") do io
        write(io, """{"zarr_format":3,"node_type":"group","attributes":{"spatialdata_attrs":{"version":"0.2.0"}}}""")
    end
    # Native metadata distinguishes this layout from Python SpatialData stores.
    open(joinpath(path, "spatialomics_meta.json"), "w") do io
        JSON.print(io, Dict(
            "format_version" => NATIVE_FORMAT_VERSION,
            "coord_systems" => Any[],
            "transforms" => Any[],
            "sources" => Any[],
        ))
    end
end

function _cleanup!(bs::BackingStore)
    bs.owned || return
    close.(values(bs.handles))
    empty!(bs.handles)
    rm(bs.path; recursive=true, force=true)
end

function _ensure_open(bs::BackingStore)
    bs.closed && throw(ArgumentError("the dataset is closed"))
    nothing
end

function _artifact_exists(bs::BackingStore, key::Tuple{Symbol,String})
    kind, name = key
    kind === :element && return any(
        isdir(joinpath(bs.path, group, name)) for group in ("points", "shapes", "images", "labels")
    )
    kind === :relation && return isdir(joinpath(bs.path, "relations", name))
    kind === :metadata && return isdir(joinpath(bs.path, "metadata", name))
    kind === :dataset && return isfile(joinpath(bs.path, "spatialomics_meta.json"))
    false
end

function _mark_dirty!(bs::BackingStore, key::Tuple{Symbol,String}; deleted::Bool=false)
    _ensure_open(bs)
    current = get(bs.changes, key, nothing)
    if deleted
        if current === :new
            delete!(bs.changes, key)
        elseif current !== :deleted
            bs.changes[key] = :deleted
        end
    elseif current !== :new
        bs.changes[key] = _artifact_exists(bs, key) ? :modified : :new
    end
    nothing
end

# ── Dataset ───────────────────────────────────────────────────────────────────

"""
    AcquisitionSource(name; region=nothing, instance_id=nothing, attributes=Dict())

A named acquisition unit such as a field of view, imaging tile, or tissue
section. `region` and `instance_id` may identify its footprint in a
`SpatialShapes` element. Observations record the source name independently of
their coordinates, so source selection remains distinct from geometric ROI
selection in overlapping acquisitions. `attributes` retains structured
technology-specific identity needed for lossless export, such as a vendor FOV
number.

# See also
[`sources`](@ref), [`source`](@ref), [`SpatialDatasetView`](@ref)
"""
struct AcquisitionSource
    name           :: String
    region_element :: Union{Nothing,String}
    region_id      :: Union{Nothing,Int32}
    attributes     :: Dict{String,Any}
end

function AcquisitionSource(name::AbstractString;
                           region::Union{Nothing,AbstractString}=nothing,
                           instance_id::Union{Nothing,Integer}=nothing,
                           attributes::AbstractDict=Dict{String,Any}())
    (region === nothing) == (instance_id === nothing) || throw(ArgumentError(
        "region and instance_id must either both be supplied or both be omitted",
    ))
    AcquisitionSource(
        String(name),
        region === nothing ? nothing : String(region),
        instance_id === nothing ? nothing : Int32(instance_id),
        deepcopy(Dict{String,Any}(
            string(key) => value for (key, value) in pairs(attributes)
        )),
    )
end

# ── Backed metadata dict ──────────────────────────────────────────────────────

"""
    BackedMetadata

Dict-like container that tracks metadata changes in its parent dataset.
Accessed as `ds.metadata`; changes become durable when [`save!`](@ref) is
called.

String and integer keys are converted to `String` automatically.
"""
mutable struct BackedMetadata <: AbstractDict{String, Any}
    data    :: Dict{String, Any}
    backing :: BackingStore
end

function Base.setindex!(bm::BackedMetadata, val, key::String)
    bm.data[key] = val
    _mark_dirty!(bm.backing, (:metadata, key))
    bm
end
Base.setindex!(bm::BackedMetadata, val, key::Union{AbstractString,Integer}) =
    setindex!(bm, val, string(key))

function Base.delete!(bm::BackedMetadata, key::String)
    haskey(bm.data, key) || throw(KeyError(key))
    delete!(bm.data, key)
    _mark_dirty!(bm.backing, (:metadata, key); deleted=true)
    bm
end
Base.delete!(bm::BackedMetadata, key::Union{AbstractString,Integer}) = delete!(bm, string(key))

Base.getindex(bm::BackedMetadata, key::String) = bm.data[key]
Base.getindex(bm::BackedMetadata, key::Union{AbstractString,Integer}) = bm[string(key)]
Base.haskey(bm::BackedMetadata, key::Union{AbstractString,Integer}) = haskey(bm.data, string(key))
Base.iterate(bm::BackedMetadata)               = iterate(bm.data)
Base.iterate(bm::BackedMetadata, state)        = iterate(bm.data, state)
Base.length(bm::BackedMetadata)                = length(bm.data)

# ── Dataset ───────────────────────────────────────────────────────────────────

"""
    SpatialDataset(; path=nothing, metadata=Dict())

Root container for a spatial omics experiment.

Holds named collections of spatial elements (`SpatialPoints`, `SpatialShapes`,
`SpatialImage`, `SpatialLabels`), a graph of `CoordinateSystem` nodes connected
by `AbstractTransformation` edges, named `SpatialRelation` objects, and free-form
metadata. Its native Zarr layout is versioned by SpatialOmics; external
SpatialData stores are handled as an import boundary.

Every dataset has a `BackingStore` Zarr directory. When `path` is `nothing`, a
temporary directory is used and cleaned up automatically. Mutations are staged
in memory and reported by [`dirty`](@ref); call [`save!`](@ref) to make them
durable. Closing a dirty dataset requires an explicit save or discard.

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
    sources       :: OrderedDict{String, AcquisitionSource}
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
        OrderedDict{String,AcquisitionSource}(),
        bs,
        Dict{String,Any}(),
        BackedMetadata(Dict{String,Any}(metadata), bs),
    )
    finalizer(ds) do d
        d.backing.closed && return
        !isempty(d.backing.changes) && @warn(
            "SpatialDataset finalized with unsaved changes",
            changes=dirty(d),
        )
        _cleanup!(d.backing)
        close.(values(d.backing.handles))
        empty!(d.backing.handles)
        d.backing.closed = true
    end
    for key in keys(metadata)
        _mark_dirty!(bs, (:metadata, string(key)))
    end
    ds
end

# ── Lifecycle ─────────────────────────────────────────────────────────────────

"""
    close(ds; discard=false)

Close a dataset and release its backing resources. A dirty dataset is rejected
unless `discard=true`; call [`save!`](@ref) or [`discard!`](@ref) first when the
changes should be kept or reviewed explicitly.
"""
function Base.close(ds::SpatialDataset; discard::Bool=false)
    ds.backing.closed && return nothing
    if isdirty(ds) && !discard
        throw(ArgumentError(
            "dataset has unsaved changes; call save!(ds), discard!(ds), or close(ds; discard=true)",
        ))
    end
    _cleanup!(ds.backing)
    ds.backing.owned = false   # prevent double-free in finalizer
    close.(values(ds.backing.handles))
    empty!(ds.backing.handles)
    ds.backing.closed = true
    nothing
end

Base.isopen(ds::SpatialDataset) = !ds.backing.closed

"""
    isdirty(ds) → Bool

Return whether `ds` contains staged changes that have not been saved.
"""
isdirty(ds::SpatialDataset) = !isempty(ds.backing.changes)

"""
    dirty(ds) → Vector{NamedTuple}

Return the staged changes in `ds`. Each entry has `kind`, `name`, and `state`
fields; `state` is `:new`, `:modified`, or `:deleted`.
"""
function dirty(ds::SpatialDataset)
    [(kind=kind, name=name, state=state) for ((kind, name), state) in ds.backing.changes]
end

"""
    touch!(ds, name) → ds

Mark a named element or relation as modified after mutation through an external
API. Prefer [`edit!`](@ref) for scoped mutation.
"""
function touch!(ds::SpatialDataset, name::String)
    if haskey(ds.elements, name)
        _mark_dirty!(ds.backing, (:element, name))
    elseif haskey(ds.relations, name)
        _mark_dirty!(ds.backing, (:relation, name))
    else
        throw(KeyError(name))
    end
    ds
end

"""
    edit!(f, ds, name)

Run `f` on a named mutable element and mark it modified before `f` is called.
The change remains dirty if `f` throws, because partial mutation may already
have occurred.

```julia
edit!(ds, "transcripts") do points
    points.feature_id[1] = 2
end
save!(ds, "transcripts")
```
"""
function edit!(f::Function, ds::SpatialDataset, name::String)
    el = ds.elements[name]
    _mark_dirty!(ds.backing, (:element, name))
    f(el)
end

"""
    keep!(ds, path=ds.backing.path) → ds

Save the dataset and mark its backing store as permanent. When `path` differs
from the current backing path, a complete snapshot is written atomically and
the dataset is rebound to it. The resulting directory is not deleted when the
dataset is closed or garbage collected.

# See also
[`with_dataset`](@ref), [`save!`](@ref)
"""
function keep!(ds::SpatialDataset, path::String=ds.backing.path)
    save!(ds; path)
    ds.backing.owned = false
    ds
end

"""
    with_dataset(f; path=nothing, kw...)

Open a temporary dataset, run `f(ds)`, then discard it and clean up the backing
store.

The dataset is always discarded in a `finally` block. Call [`keep!`](@ref) or
[`save!`](@ref) with a permanent path inside `f` when results should survive.

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
        close(ds; discard=true)
    end
end

# ── Coordinate systems + transform API ────────────────────────────────────────

function Base.push!(ds::SpatialDataset, cs::CoordinateSystem)
    ds.coord_systems[cs.name] = cs
    _mark_dirty!(ds.backing, (:dataset, "coordinate_systems"))
    ds
end

function Base.push!(ds::SpatialDataset, t::AbstractTransformation)
    push!(ds.transforms, t)
    _mark_dirty!(ds.backing, (:dataset, "coordinate_systems"))
    ds
end

function Base.push!(ds::SpatialDataset, acquisition::AcquisitionSource)
    ds.sources[acquisition.name] = acquisition
    _mark_dirty!(ds.backing, (:dataset, "coordinate_systems"))
    ds
end

"""
    elements(ds) → OrderedDict{String, Any}

Return a shallow snapshot of the named spatial elements in `ds`.

Values are concrete element types (`SpatialPoints`, `SpatialShapes`,
`SpatialImage`, `SpatialLabels`). Use the typed accessors [`points`](@ref),
[`shapes`](@ref), [`images`](@ref), [`labels`](@ref) to retrieve a specific
element with type checking.
"""
elements(ds::SpatialDataset)      = copy(ds.elements)

"""
    coord_systems(ds) → Vector{String}

Return the names of all coordinate systems registered in `ds`.
"""
coord_systems(ds::SpatialDataset) = collect(keys(ds.coord_systems))

"""
    sources(ds) → Vector{String}

Return the registered acquisition-source names in `ds`.
"""
sources(ds::SpatialDataset) = collect(keys(ds.sources))

"""
    source(ds, name) → AcquisitionSource

Return the named acquisition source. Source names can also be passed directly
to `view(ds, name)`.
"""
function source(ds::SpatialDataset, name::AbstractString)
    key = String(name)
    haskey(ds.sources, key) || throw(ArgumentError(
        "unknown acquisition source $(repr(name)); available: $(join(keys(ds.sources), ", "))",
    ))
    ds.sources[key]
end

"""
    source_attributes(source)
    source_attributes(ds, name)

Return a copy of the structured technology-specific attributes registered for
an acquisition source.
"""
source_attributes(acquisition::AcquisitionSource) = deepcopy(acquisition.attributes)
source_attributes(ds::SpatialDataset, name::AbstractString) =
    source_attributes(source(ds, name))

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
    ds.relations[name] = rel
    _mark_dirty!(ds.backing, (:relation, name))
    ds
end

function _attach_element!(ds::SpatialDataset, el, name::String)
    existing = _owning_dataset(el)
    if existing !== nothing && existing !== ds
        error("Element already attached to a different dataset. " *
              "Use `ds[\"$name\"] = copy(el)` to attach a detached copy. " *
              "Note: relations involving this element in the original dataset will not transfer.")
    end
    attachment = _dataset_ref(el)
    if existing === ds && attachment[2] != name
        error("Element is already attached to this dataset as \"$(attachment[2])\". " *
              "Use copy(el) to attach it under another name.")
    end
    if haskey(ds.elements, name) && ds.elements[name] !== el
        _clear_backref!(ds.elements[name])
    end
    _set_backref!(el, ds, name)
    ds.elements[name] = el
    _mark_dirty!(ds.backing, (:element, name))
    ds
end

function Base.delete!(ds::SpatialDataset, name::String)
    if haskey(ds.elements, name)
        el = pop!(ds.elements, name)
        _clear_backref!(el)
        _mark_dirty!(ds.backing, (:element, name); deleted=true)
    elseif haskey(ds.relations, name)
        pop!(ds.relations, name)
        _mark_dirty!(ds.backing, (:relation, name); deleted=true)
    else
        throw(KeyError(name))
    end
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
    copy(ds.relations)
end

function relations(ds::SpatialDataset, name::String)
    haskey(ds.relations, name) ||
        error("No relation \"$name\". Available: $(join(keys(ds.relations), ", "))")
    ds.relations[name]
end
