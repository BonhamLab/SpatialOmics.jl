
"""
    AbstractROI

A 2D polygon identifying a region of interest within a [`Slide`](@ref).
[`ROI`](@ref) is a generic implementation,
but specialized ROI types (eg [`FOV`](@ref))
may have additional functionality and some special properties.

Required Fields:

- loc: Geometry (eg Point or Polygon) with coordinates in slide-space

"""
abstract type AbstractROI end
abstract type AbstractExperiment end
abstract type AbstractSlide end

"""
    loc(roi::AbstractROI)

Return the location of the ROI in slide-space.

"""
loc(roi::AbstractROI) = roi.loc

"""
    struct ROI <: AbstractROI
        loc::Polygon{2} 
    end

"""
struct ROI <: AbstractROI
    loc::Polygon{2} 
end

"""
    mutable struct FOV <: AbstractROI
        roi::ROI
        img::DimArray
    end

"""
mutable struct FOV <: AbstractROI
    roi::ROI
    img::Union{Nothing, DimArray}
end

FOV(roi::ROI) = FOV(roi, nothing)

loc(fov::FOV) = loc(fov.roi)

image(fov::FOV) = fov.img
image!(fov::FOV, img) = (fov.img = img)
size(fov::FOV) = size(fov.img)
size(fov::FOV, args...) = size(fov.img, args...)

"""
    struct Slide
        parent::Experiment
        properties::LittleDict{String,Any}
        fovs::LittleDict{String,FOV}
        rois::LittleDict{String,ROI}
    end

Basic struct (type) containing information
about an individual slide (eg one run through the machine).
"""
struct Slide <: AbstractSlide
    name::String
    parent::AbstractExperiment
    props::DiskStore
    fovs::LittleDict{String,FOV}
    rois::LittleDict{String,ROI}

    function Slide(name, parent)
        slidepath = joinpath(parent.base_path, "slides", name)
        props = _build_or_get_props(slidepath, name)
        fovs = LittleDict{String,FOV}()
        rois = LittleDict{String,ROI}()
        return new(name, parent, props, fovs, rois)
    end
end

fov(slide::Slide, name) = getindex(slide.fovs, name)
fovs(slide::Slide) = keys(slide.fovs)

function fov!(slide::Slide, fov::FOV, name::String)
    setindex!(slide.fovs, fov, name)
end
fov!(slide::Slide, fov::FOV) = setindex!(slide.fovs, fov, "FOV$(lpad(length(fovs(slide)) + 1, 4, '0'))")

roi(slide::Slide, name) = getindex(slide.rois, name)
rois(slide::Slide) = keys(slide.rois)
roi!(slide::Slide, name::String, roi::ROI) = setindex!(slide.rois, name, roi)
roi!(slide::Slide, roi) = setindex!(slide.rois, "ROI$(lpad(length(rois(slide)) + 1, 3, '0'))", roi)


property(slide::Slide, name) = getindex(slide.props, name)
properties(slide::Slide) = keys(slide.props) 

function Base.display(slide::Slide)
    nfovs = length(fovs(slide))
    nrois = length(rois(slide))
    println("Slide $(slide.name) with $nfovs FOVs and $nrois ROIs")
end

"""
    Experiment(exname::String, base_path::String)

Basic struct (type) containing information about an experiment.

## Indexing

Use string indexes to pull out individual [`Slide`](@ref)s.

"""
mutable struct Experiment <: AbstractExperiment
    name::String
    base_path::Union{Nothing,String}
    slides::LittleDict{String,Slide}
    props::DiskStore

    function Experiment(name::String, base_path::String)
        props = _build_or_get_props(base_path, name)
        slides = Dict{String,Slide}()
        return new(name, base_path, slides, props)

    end
end


function Base.display(ex::Experiment)
    nslides = length(ex.slides)
    println("Experiment with $nslides slides:")
    i = 0
    for slide in values(ex.slides)
        i += 1
        println("    ", slide.name)
        if i > 5
            println("...")
        end
    end
end

Base.getindex(ex::Experiment, i) = ex.slides[i]
Base.setindex!(ex::Experiment, i, thing) = setindex!(ex.slides, i, thing)

# Iterate on experiments drops down to slides
Base.iterate(ex::Experiment) = iterate(ex.slides)
Base.iterate(ex::Experiment, state) = iterate(ex.slides, state)
Base.IteratorSize(ex::Experiment) = IteratorSize(ex.slides)
Base.length(ex::Experiment) = length(ex.slides)
Base.isdone(ex::Experiment) = isdone(ex.slides)
Base.isdone(ex::Experiment, state) = isdone(ex.slides, state)
