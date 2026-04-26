# ── SpatialImage ──────────────────────────────────────────────────────────────

mutable struct SpatialImage{T, N}
    data          :: AbstractArray{T, N}
    pyramid       :: Vector{AbstractArray{T, N}}   # coarser levels; empty until built
    axes          :: NTuple{N, Symbol}              # e.g. (:c,:y,:x) or (:y,:x)
    channel_names :: Vector{String}
    coord_system  :: String
    pixel_to_cs   :: AbstractTransformation        # image pixel coords → coord_system
end

function _default_image_axes(N::Int)
    N == 2 ? (:y, :x) :
    N == 3 ? (:c, :y, :x) :
    ntuple(i -> Symbol("dim_$i"), N)
end

function SpatialImage(data::AbstractArray{T, N};
                      axes          = _default_image_axes(N),
                      channel_names = String[],
                      coord_system  = "",
                      pixel_to_cs   = Identity("pixel", coord_system)) where {T, N}
    SpatialImage{T, N}(
        data,
        AbstractArray{T, N}[],
        NTuple{N, Symbol}(axes),
        channel_names,
        coord_system,
        pixel_to_cs)
end

# ── Accessors ──────────────────────────────────────────────────────────────────

data(img::SpatialImage)          = img.data
channel_names(img::SpatialImage) = img.channel_names
coord_system(img::SpatialImage)  = img.coord_system

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

function build_pyramid!(img::SpatialImage, n_levels::Int=3)
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

struct SpatialLabels{T<:Integer, N}
    data         :: AbstractArray{T, N}
    axes         :: NTuple{N, Symbol}
    instance_map :: Dict{T, Int32}       # pixel label → canonical instance_id
    coord_system :: String
    pixel_to_cs  :: AbstractTransformation
end

function SpatialLabels(data::AbstractArray{T, N};
                       axes         = _default_image_axes(N),
                       instance_map = Dict{T, Int32}(),
                       coord_system = "",
                       pixel_to_cs  = Identity("pixel", coord_system)) where {T<:Integer, N}
    SpatialLabels{T, N}(data, NTuple{N, Symbol}(axes), instance_map, coord_system, pixel_to_cs)
end

coord_system(lbl::SpatialLabels)  = lbl.coord_system
instance_ids(lbl::SpatialLabels)  = sort(unique(values(lbl.instance_map)))
Base.size(lbl::SpatialLabels)     = size(lbl.data)

# ── Dataset typed accessors ────────────────────────────────────────────────────

function images(ds::SpatialDataset, name::String)
    el = ds.elements[name]
    el isa SpatialImage || error("Element \"$name\" is not SpatialImage (got $(typeof(el)))")
    el
end

function labels(ds::SpatialDataset, name::String)
    el = ds.elements[name]
    el isa SpatialLabels || error("Element \"$name\" is not SpatialLabels (got $(typeof(el)))")
    el
end
