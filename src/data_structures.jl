

"""
    AbstractROI

A 2D polygon identifying a region of interest within a [`Slide`](@ref).
[`ROI`](@ref) is a generic implementation,
but specialized ROI types (eg [`FOV`](@ref))
may have additional functionality and some special properties.

Required Fields:

- name: A string used to identify the ROI.
  All ROI names within a slide must have unique names
- parent: [`Slide`](@ref) that contains the ROI.
- loc: Geometry (eg Point or Polygon) with coordinates in slide-space
- props: Disk-written metadata dictionary for storing arbitrary information

"""
abstract type AbstractROI end
abstract type AbstractExperiment end
abstract type AbstractSlide end

struct ROI <: AbstractROI
    name::String
    parent::AbstractSlide
    loc::Polygon{2} 
    props::DiskStore
end

"""
    mutable struct FOV <: AbstractROI
        parent::Slide
        loc::Rectf{2} 
        props::LittleDict{String,Any}
    end

"""
mutable struct FOV <: AbstractROI
    parent::AbstractSlide
    loc::Rectf{2} 
    props::DiskStore
end


function Base.display(fov::FOV)
    sn = string(fov.parent.name)
    println("FOV $(fov.name)  in slide $sn")
end

"""
    struct Slide
        parent::Experiment
        properties::LittleDict{String,Any}
        rois::LittleDict{String,AbstractROI}
    end

Basic struct (type) containing information
about an individual slide (eg one run through the machine).

"""
struct Slide <: AbstractSlide
    parent::AbstractExperiment
    props::DiskStore
    fovs::LittleDict{String,FOV}
    rois::LittleDict{String,ROI}
    
    function Slide(name, parent)
        slidepath = joinpath(parent.base_path, name)
        props = _build_or_get_props(slidepath, name)
        fovs = LittleDict{String,FOV}()
        rois = LittleDict{String,ROI}()
        return new(parent, props, fovs, rois)
    end
end

roi(slide::Slide, name) = getindex(slide.rois, name)
rois(slide::Slide) = keys(slide.rois) 

property(slide::Slide, name) = getindex(slide.props, name)
properties(slide::Slide) = keys(slide.props) 

function Base.display(slide::Slide)
    nrois = sum(g -> length(rois(slide)))
    println("Slide $(slide.name) with $nrois ROIs")
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
        return new(name, base_path, slides)

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


