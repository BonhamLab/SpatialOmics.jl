# ── SpatialImage ──────────────────────────────────────────────────────────────

"""
    SpatialImage{T, N}

N-dimensional raster image in a named coordinate system, optionally with
pre-computed pyramid levels.

The image data array has axes labelled by `axes` — a tuple of symbols such as
`(:c, :y, :x)` for a 3-D multi-channel image or `(:y, :x)` for a 2-D image.
The `pixel_to_cs` transformation maps pixel-space indices to the physical
coordinate system. Pyramid levels (coarser sub-sampled copies) are stored in
`pyramid` and can be built with `build_pyramid!`.

A `display_transform` (e.g. from `scaleminmax`) is applied after data
materialisation at render time, never wrapped around Zarr arrays.

# Constructors
    SpatialImage(data; axes, channel_names, coord_system, pixel_to_cs, pyramid, display_transform)

# See also
[`SpatialLabels`](@ref), [`channel`](@ref), [`scaleminmax`](@ref),
[`build_pyramid!`](@ref), [`colorview`](@ref)
"""
mutable struct SpatialImage{T, N}
    data              :: AbstractArray{T, N}
    pyramid           :: Vector{AbstractArray{T, N}}   # coarser levels; empty until built
    axes              :: NTuple{N, Symbol}              # e.g. (:c,:y,:x) or (:y,:x)
    channel_names     :: Vector{String}
    coord_system      :: String
    pixel_to_cs       :: AbstractTransformation        # image pixel coords → coord_system
    display_transform :: Union{Nothing, Function}      # applied post-materialization in display
    _attachment       :: Union{Nothing, Tuple{WeakRef, String}}
end

function _default_image_axes(N::Int)
    N == 2 ? (:y, :x) :
    N == 3 ? (:c, :y, :x) :
    ntuple(i -> Symbol("dim_$i"), N)
end

function SpatialImage(data::AbstractArray{T, N};
                      axes              = _default_image_axes(N),
                      channel_names     = String[],
                      coord_system      = "",
                      pixel_to_cs       = Identity("pixel", coord_system),
                      pyramid           = nothing,
                      display_transform = nothing) where {T, N}
    pyr = pyramid === nothing ? AbstractArray{T,N}[] :
          AbstractArray{T,N}[p for p in pyramid]
    SpatialImage{T, N}(data, pyr, NTuple{N, Symbol}(axes), channel_names,
                       coord_system, pixel_to_cs, display_transform, nothing)
end

_dataset_ref(img::SpatialImage) = img._attachment
function _set_backref!(img::SpatialImage, ds::SpatialDataset, name::String)
    img._attachment = (WeakRef(ds), name)
end
function _clear_backref!(img::SpatialImage)
    img._attachment = nothing
    img
end

Base.setindex!(ds::SpatialDataset, img::SpatialImage, name::String) =
    _attach_element!(ds, img, name)

# ── Accessors ──────────────────────────────────────────────────────────────────

"""
    data(img) → AbstractArray

Return the underlying data array of a `SpatialImage` or `SpatialLabels`.
May be a Zarr-backed `DiskArray` — use `Array(data(img))` to force a full read.
"""
data(img::SpatialImage)          = img.data

"""
    channel_names(img) → Vector{String}

Return the channel names of a `SpatialImage`. Empty if names were not set.

# See also
[`nchannels`](@ref), [`channel`](@ref)
"""
channel_names(img::SpatialImage) = img.channel_names
coord_system(img::SpatialImage)  = img.coord_system

"""
    nchannels(img) → Int

Return the number of channels in a `SpatialImage`. Returns `1` for 2-D images
with no channel axis.

# See also
[`channel_names`](@ref), [`channel`](@ref)
"""
function nchannels(img::SpatialImage)
    idx = findfirst(==(:c), img.axes)
    idx === nothing ? 1 : size(img.data, idx)
end

Base.size(img::SpatialImage)   = size(img.data)
Base.length(img::SpatialImage) = length(img.data)

# ── Pyramid ────────────────────────────────────────────────────────────────────

function _spatial_dims(axes::NTuple{N, Symbol}) where N
    Tuple(i for (i, a) in enumerate(axes) if a in (:x, :y, :z))
end

"""
    build_pyramid!(img, n_levels=3) → img

Build a Gaussian downsampling pyramid in-place, storing `n_levels` progressively
coarser arrays in `img.pyramid`.

Each level halves the spatial resolution along the `:x` and `:y` axes using
`ImageBase.restrict`. The channel axis (`:c`) is not downsampled. Existing
pyramid levels are discarded before building.

# See also
[`scaleminmax`](@ref), [`channel`](@ref)
"""
function build_pyramid!(img::SpatialImage, n_levels::Int=3)
    owner = _owning_dataset(img)
    owner === nothing || touch!(owner, _dataset_ref(img)[2])
    empty!(img.pyramid)
    sdims   = _spatial_dims(img.axes)
    current = img.data
    for _ in 1:n_levels
        current = restrict(current, sdims)
        push!(img.pyramid, current)
    end
    img
end

# ── SpatialLabels ─────────────────────────────────────────────────────────────

"""
    SpatialLabels{T<:Integer, N}

Integer segmentation mask in a named coordinate system, where each pixel value
identifies an object instance.

The `instance_map` dictionary maps raw pixel label values to canonical
`instance_id` values, matching the convention used in `SpatialPoints` and
`SpatialShapes`. Use `instance_ids(lbl)` to list all non-background instances.

# Constructors
    SpatialLabels(data; axes, instance_map, coord_system, pixel_to_cs)

# See also
[`SpatialImage`](@ref), [`instance_ids`](@ref), [`data`](@ref)
"""
struct SpatialLabels{T<:Integer, N}
    data         :: AbstractArray{T, N}
    axes         :: NTuple{N, Symbol}
    instance_map :: Dict{T, Int32}       # pixel label → canonical instance_id
    coord_system :: String
    pixel_to_cs  :: AbstractTransformation
end

Base.setindex!(ds::SpatialDataset, lbl::SpatialLabels, name::String) =
    _attach_element!(ds, lbl, name)

function SpatialLabels(data::AbstractArray{T, N};
                       axes         = _default_image_axes(N),
                       instance_map = Dict{T, Int32}(),
                       coord_system = "",
                       pixel_to_cs  = Identity("pixel", coord_system)) where {T<:Integer, N}
    SpatialLabels{T, N}(data, NTuple{N, Symbol}(axes), instance_map, coord_system, pixel_to_cs)
end

coord_system(lbl::SpatialLabels)  = lbl.coord_system
data(lbl::SpatialLabels)          = lbl.data

"""
    instance_ids(lbl) → Vector{Int32}

Return the sorted list of unique instance IDs present in a `SpatialLabels` mask.

Does not include the background (pixels not in `instance_map`).

# See also
[`instance_id`](@ref), [`SpatialLabels`](@ref)
"""
instance_ids(lbl::SpatialLabels)  = sort(unique(values(lbl.instance_map)))
Base.size(lbl::SpatialLabels)     = size(lbl.data)

# ── Dataset typed accessors ────────────────────────────────────────────────────

"""
    images(ds, name) → SpatialImage
    images(v, name) → SpatialImage  (cropped to view extent)

Retrieve the named `SpatialImage` from a dataset or dataset view.
Raises an error if the element is not a `SpatialImage`.

# See also
[`labels`](@ref), [`points`](@ref), [`shapes`](@ref)
"""
function images(ds::SpatialDataset, name::String)
    el = ds.elements[name]
    el isa SpatialImage || error("Element \"$name\" is not SpatialImage (got $(typeof(el)))")
    el
end

"""
    labels(ds, name) → SpatialLabels

Retrieve the named `SpatialLabels` from a dataset or dataset view.
Raises an error if the element is not a `SpatialLabels`.

# See also
[`images`](@ref), [`instance_ids`](@ref)
"""
function labels(ds::SpatialDataset, name::String)
    el = ds.elements[name]
    el isa SpatialLabels || error("Element \"$name\" is not SpatialLabels (got $(typeof(el)))")
    el
end

# ── channel — lazy 2D slice extraction ────────────────────────────────────────

"""
    channel(img, ch) → SpatialImage

Return a lazy view of a single channel from a multi-channel `SpatialImage`.

`ch` can be an `Int` (1-based channel index) or a `String` (channel name from
`channel_names(img)`). The result is a 2-D `SpatialImage` with the channel axis
removed. Pyramid levels are sliced correspondingly.

# See also
[`nchannels`](@ref), [`channel_names`](@ref), [`colorview`](@ref)
"""
function channel(img::SpatialImage{T,N}, ch::Int) where {T,N}
    ci = findfirst(==(:c), img.axes)
    ci === nothing && return img
    slices   = ntuple(d -> d == ci ? ch : Colon(), N)
    new_axes = NTuple{N-1, Symbol}(img.axes[d] for d in 1:N if d != ci)
    SpatialImage(view(img.data, slices...);
                 axes=new_axes, coord_system=img.coord_system,
                 pixel_to_cs=img.pixel_to_cs,
                 pyramid=AbstractArray{T, N-1}[view(lvl, slices...) for lvl in img.pyramid],
                 display_transform=img.display_transform)
end

function channel(img::SpatialImage, ch::String)
    i = findfirst(==(ch), img.channel_names)
    i === nothing && error("Channel \"$ch\" not found. Available: $(img.channel_names)")
    channel(img, i)
end

# ── scaleminmax — lazy display-time intensity rescaling ────────────────────────

"""
    scaleminmax(img) → SpatialImage

Return a copy of `img` with a min-max intensity rescaling transform set as its
`display_transform`.

The intensity range is sampled from the coarsest pyramid level (a fast single
bulk read). The transform is applied after Zarr materialisation at display
time — it is never wrapped around lazy disk arrays.

# See also
[`channel`](@ref), [`colorview`](@ref), [`build_pyramid!`](@ref)
"""
function scaleminmax(img::SpatialImage)
    src      = isempty(img.pyramid) ? img.data : img.pyramid[end]
    mn, mx   = Float32.(extrema(Array(src)))
    SpatialImage(img.data;
                 axes=img.axes, channel_names=img.channel_names,
                 coord_system=img.coord_system, pixel_to_cs=img.pixel_to_cs,
                 pyramid=img.pyramid,
                 display_transform=scaleminmax(mn, mx))
end

# ── pyramid_level — internal helper (not exported) ─────────────────────────────

function _pyramid_level(img::SpatialImage, level::Int)
    level == 0 && return img
    1 <= level <= length(img.pyramid) ||
        error("Level $level out of range [0, $(length(img.pyramid))]")
    SpatialImage(img.pyramid[level];
                 axes          = img.axes,
                 channel_names = img.channel_names,
                 coord_system  = img.coord_system,
                 pixel_to_cs   = _scale_pixel_transform(img.pixel_to_cs, 2.0^level),
                 display_transform = img.display_transform)
end

_scale_pixel_transform(t::Identity, scale::Float64) =
    Affine(SMatrix{3,3,Float64}(scale,0,0, 0,scale,0, 0,0,1), t.src, t.dst)

function _scale_pixel_transform(t::Affine, scale::Float64)
    S = SMatrix{3,3,Float64}(scale,0,0, 0,scale,0, 0,0,1)
    Affine(t.matrix * S, t.src, t.dst)
end

_scale_pixel_transform(t::Sequence, scale::Float64) =
    Sequence([_scale_pixel_transform(t.steps[1], scale); t.steps[2:end]], t.src, t.dst)

# ── SpatialImageColorView ─────────────────────────────────────────────────────

"""
    SpatialImageColorView{C<:Colorant, T, N}

Display-only wrapper that pairs a `SpatialImage` with a colorant type and an
optional scalar display transform.

Produced by `colorview(CT, img)`, `colorview(CT, img, ch)`, or
`colorview(CT, img1, img2, ...)` for RGB composites. The `display_transform`
(if any) is applied after Zarr materialisation, never wrapping disk arrays.

# See also
[`colorview`](@ref), [`scaleminmax`](@ref), [`channel`](@ref)
"""
struct SpatialImageColorView{C<:Colorant, T, N}
    data         :: AbstractArray{T, N}          # raw data (may be zarr-backed)
    pyramid      :: Vector{AbstractArray{T, N}}  # raw pyramid levels (may be zarr-backed)
    colorant     :: Type{C}                      # colorant applied after materialisation
    transform    :: Union{Nothing, Function}     # scalar map applied before colorview
    coord_system :: String
    pixel_to_cs  :: AbstractTransformation
    axes         :: NTuple{N, Symbol}
end

# ── colorview extensions ───────────────────────────────────────────────────────

"""
    colorview(CT, img) → SpatialImageColorView
    colorview(CT, img, ch) → SpatialImageColorView
    colorview(CT, img1, img2, ...) → SpatialImageColorView

Wrap one or more `SpatialImage` objects with colorant type `CT` for display.

Single-image forms extract `channel(img, 1)` if `img` is multi-channel unless
`ch` is specified. The multi-image form creates an RGB (or N-channel) composite;
all inputs must be 2-D (call `channel` first).

```julia
cview = colorview(Gray, scaleminmax(channel(img, 1)))
```

# See also
[`SpatialImageColorView`](@ref), [`channel`](@ref), [`scaleminmax`](@ref)
"""
function colorview(CT::Type{<:Colorant}, img::SpatialImage)
    ndims(img.data) > 2 ? _spatial_colorview(CT, channel(img, 1)) :
                          _spatial_colorview(CT, img)
end

colorview(CT::Type{<:Colorant}, img::SpatialImage, ch) =
    _spatial_colorview(CT, channel(img, ch))

function colorview(CT::Type{<:Colorant}, imgs::SpatialImage...)
    all(si -> ndims(si.data) == 2, imgs) ||
        error("All SpatialImage arguments must be 2D for composite colorview; call channel() first")
    _tf(si, arr) = isnothing(si.display_transform) ? arr : mappedarray(si.display_transform, arr)
    data_arrs = [_tf(si, si.data) for si in imgs]
    cdata     = colorview(CT, data_arrs...)          # ImageCore dispatch — lazy colored composite
    C         = eltype(cdata)
    n_levels  = isempty(imgs[1].pyramid) ? 0 : minimum(length(si.pyramid) for si in imgs)
    cpyr      = AbstractArray{C, 2}[
        colorview(CT, [_tf(imgs[j], imgs[j].pyramid[i]) for j in 1:length(imgs)]...) for i in 1:n_levels
    ]
    SpatialImageColorView{C, C, 2}(cdata, cpyr, C, nothing,
                                   imgs[1].coord_system, imgs[1].pixel_to_cs, imgs[1].axes)
end

function _spatial_colorview(CT::Type{<:Colorant}, img::SpatialImage{T,N}) where {T,N}
    SpatialImageColorView{CT, T, N}(img.data, img.pyramid, CT, img.display_transform,
                                    img.coord_system, img.pixel_to_cs, img.axes)
end

# ── Base.view — lazy spatial crop of SpatialImage to a SpatialExtent ──────────
# Maps extent corners through inverse of pixel_to_cs to get pixel indices,
# returns a SpatialImage with view-sliced data and adjusted pixel_to_cs.

function Base.view(img::SpatialImage, ext::SpatialExtent)
    N  = ndims(img.data)
    xi = something(findfirst(==(:x), img.axes), 1)
    yi = something(findfirst(==(:y), img.axes), 2)
    lo = _global_to_pixel(img.pixel_to_cs, SVector(ext.xmin, ext.ymin))
    hi = _global_to_pixel(img.pixel_to_cs, SVector(ext.xmax, ext.ymax))
    xi_lo, xi_hi = _px_range(lo[1], hi[1], size(img.data, xi))
    yi_lo, yi_hi = _px_range(lo[2], hi[2], size(img.data, yi))
    sl   = ntuple(d -> d == xi ? (xi_lo:xi_hi) : d == yi ? (yi_lo:yi_hi) : Colon(), N)
    p2cs = _shift_pixel_origin(img.pixel_to_cs, Float64(xi_lo-1), Float64(yi_lo-1))
    nx, ny = size(img.data, xi), size(img.data, yi)
    new_pyr = map(img.pyramid) do lvl
        sx = nx / size(lvl, xi)   # actual scale — works for build_pyramid! and pre-stored Zarr levels
        sy = ny / size(lvl, yi)
        pl, ph = _px_range(lo[1]/sx, hi[1]/sx, size(lvl, xi))
        ql, qh = _px_range(lo[2]/sy, hi[2]/sy, size(lvl, yi))
        view(lvl, ntuple(d -> d == xi ? (pl:ph) : d == yi ? (ql:qh) : Colon(), N)...)
    end
    SpatialImage(view(img.data, sl...);
                 axes=img.axes, channel_names=img.channel_names,
                 coord_system=img.coord_system, pixel_to_cs=p2cs,
                 pyramid=new_pyr, display_transform=img.display_transform)
end

function Base.view(lbl::SpatialLabels, ext::SpatialExtent)
    N  = ndims(lbl.data)
    xi = something(findfirst(==(:x), lbl.axes), 1)
    yi = something(findfirst(==(:y), lbl.axes), 2)
    lo = _global_to_pixel(lbl.pixel_to_cs, SVector(ext.xmin, ext.ymin))
    hi = _global_to_pixel(lbl.pixel_to_cs, SVector(ext.xmax, ext.ymax))
    xi_lo, xi_hi = _px_range(lo[1], hi[1], size(lbl.data, xi))
    yi_lo, yi_hi = _px_range(lo[2], hi[2], size(lbl.data, yi))
    slices = ntuple(
        dimension -> dimension == xi ? (xi_lo:xi_hi) :
                     dimension == yi ? (yi_lo:yi_hi) : Colon(),
        N,
    )
    pixel_to_cs = _shift_pixel_origin(
        lbl.pixel_to_cs, Float64(xi_lo - 1), Float64(yi_lo - 1),
    )
    SpatialLabels(
        view(lbl.data, slices...);
        axes=lbl.axes,
        instance_map=copy(lbl.instance_map),
        coord_system=lbl.coord_system,
        pixel_to_cs,
    )
end

_px_range(lo, hi, n) = (clamp(floor(Int, min(lo, hi)) + 1, 1, n),
                         clamp(ceil(Int,  max(lo, hi)),     1, n))

_global_to_pixel(::Identity, p::SVector) = p
_global_to_pixel(t::Affine,  p::SVector) =
    SVector{2,Float64}((SMatrix{3,3,Float64}(inv(Matrix(t.matrix))) *
                         SVector(p[1], p[2], 1.0))[1:2])

_shift_pixel_origin(t::Identity, dx::Float64, dy::Float64) =
    (iszero(dx) && iszero(dy)) ? t : translation(dx, dy, t.src, t.dst)

function _shift_pixel_origin(t::Affine, dx::Float64, dy::Float64)
    (iszero(dx) && iszero(dy)) && return t
    S = SMatrix{3,3,Float64}(1,0,0, 0,1,0, dx,dy,1)
    Affine(t.matrix * S, t.src, t.dst)
end
