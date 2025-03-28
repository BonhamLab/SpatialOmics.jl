
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


"""
    struct Slide
        name::String
        parent::Experiment
        props::LittleDict{String,Any}
        fovs::LittleDict{String,FOV}
        rois::LittleDict{String,ROI}
    end

Basic struct (type) containing information
about an individual slide (eg one run through the machine).
"""
struct Slide <: AbstractSlide
    name::String
    parent::AbstractExperiment
    props::LittleDict{String,Any}
    fovs::LittleDict{String,FOV}
    rois::LittleDict{String,ROI}
end



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
    props::LittleDict{String,Any}
    Experiment(name::String) = new(name, nothing, LittleDict{String,Slide}(), LittleDict{String,Any}())
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

## Path handling

experiment_path(ex::Experiment) = joinpath(ex.base_path, ex.name)
experiment_path(slide::Slide) = experiment_path(slide.parent)

slides_path(ex::Experiment) = joinpath(experiment_path(ex), "slides")
slides_path(slide::Slide) = slides_path(slide.parent)

slide_path(ex::Experiment, slidename::String) = joinpath(slides_path(ex), slidename)
slide_path(slide::Slide) = joinpath(slides_path(slide), slide.name)

rois_path(ex::Experiment, slidename::String) = joinpath(slide_path(ex, slidename), "rois")
fovs_path(ex::Experiment, slidename::String) = joinpath(slide_path(ex, slidename), "fovs")
imgs_path(ex::Experiment, slidename::String) = joinpath(slide_path(ex, slidename), "imgs")

rois_path(slide::Slide) = joinpath(slide_path(slide), "rois")
fovs_path(slide::Slide) = joinpath(slide_path(slide), "fovs")
imgs_path(slide::Slide) = joinpath(slide_path(slide), "imgs")

fov_path(slide::Slide, fovname::String) = joinpath(fovs_path(slide), "fovs", "$fovname.toml")
roi_path(slide::Slide, roiname::String) = joinpath(rois_path(slide), "rois", "$roiname.toml")
img_path(slide::Slide, imgname::String; ext="tiff") = joinpath(imgs_path(slide), "imgs", "$imgname.$ext")

properties_path(ex::Experiment) = joinpath(experiment_path, "$(ex.name)_properties.toml")
properties_path(slide::Slide) = joinpath(slide_path(slide), "$(ex.name)_properties.toml")

## Constructors

function ROI(properties::String)
    props = TOML.parsefile(properties)
    loc = Polygon(Point2f.([
        (point["x"], point["y"]) for point in props["points"]
    ]))
    return ROI(loc)
end

function FOV(properties::String)
    roi = ROI(properties)
    return FOV(roi, nothing)
end

"""
    roi!(slide::Slide[, name::String], roi::ROI)

Adds ROI to a slide, and saves the points of the ROI to disk.
The expected format of the ROI is a TOML file with a `points` key
containing a list of points, eg:

```toml
[[points]]  # First point
x = 0.0
y = 0.0

[[points]]  # Second point
x = 1.0
y = 1.0
#...
```

If `name` is not provided, the ROI is named "ROI001", "ROI002", etc.,
based on the number of ROIs already present in the slide.
"""
function roi!(slide::Slide, name::String, roi::ROI)
    rpath = roi_path(slide, name)
    open(rpath, "w") do io
        TOML.print(io, Dict(
            "points" => [
                Dict("x" => point[1], "y" => point[2]) for point in points(roi)
            ]
        ))
    end
    setindex!(slide.rois, name, roi)
end

roi!(slide::Slide, roi::ROI) = roi!(slide.rois, "ROI$(lpad(length(rois(slide)) + 1, 3, '0'))", roi)

"""
    roi!(slide::Slide, properties::String)

Loads ROI from disk and adds it to a [`Slide`](@ref).
The expected format of `properties` is a TOML file with a `points` key
containing a list of points, eg:

```toml
[[points]]  # First point
x = 0.0
y = 0.0

[[points]]  # Second point
x = 1.0
y = 1.0
#...
```

The name of the ROI is taken from the filename.
"""
function roi!(slide::Slide, properties::String)
    name = splitext(basename(properties))[1]
    roi = ROI(properties)
    setindex!(slide.rois, name, roi)
end

"""
    fov!(slide::Slide[, name::String], fov::FOV)

Adds FOV to a slide, and saves the points of the FOV to disk.
The expected format of the FOV is a TOML file with a `points` key
containing a list of points, eg: 

```toml
[[points]]  # First point
x = 0.0
y = 0.0

[[points]]  # Second point
x = 1.0
y = 1.0
#...
```

If `name` is not provided, the FOV is named "FOV001", "FOV002", etc.,
based on the number of FOVs already present in the slide.
"""
function fov!(slide::Slide, name::String, fov::FOV)
    fpath = fov_path(slide, name)
    open(fpath, "w") do io
        TOML.print(io, Dict(
            "points" => [
                Dict("x" => point[1], "y" => point[2]) for point in points(fov)
            ]
        ))
    end
    setindex!(slide.fovs, name, fov)
end

fov!(slide::Slide, fov::FOV) = fov!(slide.fovs, "FOV$(lpad(length(fovs(slide)) + 1, 4, '0'))", fov)

"""
    fov!(slide::Slide, properties::String)

Loads FOV from disk and adds it to a [`Slide`](@ref).
The expected format of `properties` is a TOML file with a `points` key
containing a list of points, eg:

```toml
[[points]]  # First point
x = 0.0
y = 0.0

[[points]]  # Second point
x = 1.0
y = 1.0
#...
```

The name of the FOV is taken from the filename.
"""
function fov!(slide::Slide, properties::String)
    name = splitext(basename(properties))[1]
    fov = FOV(properties)
    setindex!(slide.fovs, name, fov)
end


function load_experiment!(ex::Experiment, base_path::String)
    ex.base_path = base_path
    for slidedir in readdir(slides_path(ex); join=true)
        slidename = basename(slidedir)
        isfile(joinpath(slidedir, "$(slidename)_properties.toml")) || continue
        props = TOML.parsefile(joinpath(slidedir, "$(slidename)_properties.toml")) |> LittleDict
               
        slide = Slide(slidedir, ex, props, LittleDict{String,FOV}(), LittleDict{String,ROI}())
        for fov in readdir(fovs_path(slide))
            fov!(slide, joinpath(fovs_path(slide), fov))
        end
        for roi in readdir(rois_path(slide))
            roi!(slide, joinpath(rois_path(slide), roi))
        end
        setindex!(ex.slides, slide.name, slide)
    end
end

## Accessors

"""
    loc(roi::AbstractROI)

Return the location of the ROI in slide-space.

"""
loc(roi::AbstractROI) = roi.loc

points(roi::AbstractROI) = loc(roi).exterior

loc(fov::FOV) = loc(fov.roi)

image(fov::FOV) = fov.img
image!(fov::FOV, img) = (fov.img = img)
size(fov::FOV) = size(fov.img)
size(fov::FOV, args...) = size(fov.img, args...)


fov(slide::Slide, name) = getindex(slide.fovs, name)
fovs(slide::Slide) = keys(slide.fovs)

function fov!(slide::Slide, fov::FOV, name::String) 
    setindex!(slide.fovs, fov, name)
end

fov!(slide::Slide, fov::FOV) = fov!(slide, fov, "FOV$(lpad(length(fovs(slide)) + 1, 4, '0'))")

roi(slide::Slide, name) = getindex(slide.rois, name)
rois(slide::Slide) = keys(slide.rois)

property(slide::Slide, name) = getindex(slide.props, name)
properties(slide::Slide) = keys(slide.props) 
property!(slide::Slide, name, value) = setindex!(slide.props, value, name)