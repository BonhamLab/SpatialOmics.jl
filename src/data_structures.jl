abstract type AbstractSlide end
abstract type AbstractImageGroup end
abstract type AbstractFieldOfView end


"""
    Experiment(exname::String, base_path::String)

Basic struct (type) containing information about an experiment.
Expects a folder prepared by [`prepare_cosmx_experiment!`](@ref).

## Indexing

Use string indexes to pull out individual [`Slide`](@ref)s.

```julia-repl
julia> julia> ex = SO.Experiment("tlr78", "data/processed")
[2024-10-18 16:40:16] Info: ["mw_mus_p1_07", "mw_mus_p1_09"]
Experiment with 2 slides:
    mw_mus_p1_07
    mw_mus_p1_09

julia> ex["mw_mus_p1_09"]
Slide mw_mus_p1_09 with 30 imagegroups and 191 FOVs
```
"""
mutable struct Experiment{T}
    name::String
    base_path::Union{Nothing,String}
    remote
    slides::LittleDict{String,T}

    function Experiment{T}(name::String, base_path::String) where {T<:AbstractSlide}
        slides = Dict{String,T}()
        return new{T}(name, base_path, nothing, slides)

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

Base.dirname(ex::Experiment) = joinpath(ex.base_path, ex.name)


"""
    struct Slide{T}
        name::String
        properties::LittleDict{String,Any}
        experiment::Experiment
        positions
        transcripts
        imagegroups::Matrix{T}
    end

Basic struct (type) containing information
about an individual slide (eg one run through the machine).

Numerical or string indexing (#TODO) will
pull out individual image imagegroups.
using [`display_slide`](@ref) with Makie loaded
to view [imagegroups](@ref ImageGroup) and their individual [FOVs](@ref FOV).

"""
struct Slide{T} <: AbstractSlide
    name::String
    parent::Experiment
    properties::LittleDict
    positions::DataFrame
    transcripts::GroupedDataFrame
    imagegroups::LittleDict
end

Base.getindex(cs::Slide, imagegroup_i) = cs.imagegroups[imagegroup_i]

function Base.display(cs::Slide)
    nimagegroups = length(cs.imagegroups)
    nfovs = sum(g -> length(g.fovs), values(cs.imagegroups))
    println("Slide $(cs.name) with $nimagegroups imagegroups and $nfovs FOVs")
end


"""
    struct ImageGroup{T}
        id::Int16
        parent::Slide
        fovs::Matrix{T}
    end

Superstructure of [`FOV`](@ref)s,
containing FOVs captured together.
Most mutating methods applying to FOVs can also be performed on
`ImageGroup`s, and will affect all underlying FOVs.
"""
struct ImageGroup{T} <: AbstractImageGroup
    id::Int16
    parent::Slide
    fovs::Matrix{T}
end

Base.getindex(imagegroup::ImageGroup, inds...) = getindex(imagegroup.fovs, inds...)

function Base.display(imagegroup::ImageGroup)
    println("ImageGroup with FOVs:")
    display([f.id for f in imagegroup.fovs])
    if all(f -> !isnothing(f.thumbnail), imagegroup.fovs)
        display(hvcat(size(imagegroup.fovs, 2), (colorview(RGB, (f.thumbnail[:, :, i] for i in 3:5)...) for f in permutedims(imagegroup.fovs))...))
    end
end

"""
    mutable struct FOV
        id::Int16
        parent::ImageGroup
        img::Union{Nothing, <:AbstractArray}
        thumbnail::Union{Nothing, <:AbstractArray}
    end

Most basic unit of a [`Experiment`](@ref).
May contain information about underlying image and thumbnail
(see also [`prepare_image`](@ref), [`fov_image`](@ref), and [`fov_image!`](@ref)).
"""
mutable struct FOV <: AbstractFieldOfView
    id::Int16
    parent::ImageGroup
    img::Union{Nothing,<:AbstractArray}
    thumbnail::Union{Nothing,<:AbstractArray}
end

function Base.display(fov::FOV)
    println("FOV $(fov.id) from imagegroup $(fov.parent.id) in slide $(fov.parent.slide.name)")
    if !isnothing(fov.thumbnail)
        display(_cosmx_colorview(fov.thumbnail))
    elseif !isnothing(fov.img)
        display(_cosmx_colorview(fov.img))
    end
end

################
# Constructors #
################

FOV(id::Integer, parent::ImageGroup) = FOV(id, parent, nothing, nothing)

"""
    ImageGroup(imagegroup_id, slide)


"""
function ImageGroup(imagegroup_id, slide)
    imagegroup_def = slide.properties["imagegroups"][imagegroup_id]
    return ImageGroup{FOV}(imagegroup_id, slide, Matrix{FOV}(undef, size(imagegroup_def)...))
end

"""
    Slide(slname, experiment)

Constructor for [`Slide`](@ref)
"""
function Slide(slname, experiment::Experiment)
    experiment_path = dirname(experiment)
    slide_files = _slide_file_names(experiment_path, slname)

    properties = YAML.load_file(slide_files.properties)
    _normalize_imagegroups!(properties)
    positions = copy(Arrow.Table(slide_files.positions) |> DataFrame)

    transcripts = groupby(Arrow.Table(slide_files.transcripts) |> DataFrame, "fov")

    slide = Slide{ImageGroup}(slname, experiment, properties, positions, transcripts, LittleDict{Any,ImageGroup}())

    pos_imagegroups = groupby(positions, "imagegroup_id")

    for (imagegroup_id, fov_ids) in properties["imagegroups"]
        imagegroup = ImageGroup(imagegroup_id, slide)
        imagegroupdf = pos_imagegroups[(; imagegroup_id)]
        for I in eachindex(fov_ids)
            fov_id = fov_ids[I]
            row = only(subset(imagegroupdf, "fov" => ByRow(==(fov_id))))
            if ismissing(row.dir) || ismissing(row.file)
                @warn "$imagegroup_id does not have an associated image"
                fovpath = ""
            else
                fovpath = joinpath(row.dir, row.file)
            end

            fov = FOV(fov_id, imagegroup)
            imagegroup.fovs[I] = fov
        end

        slide.imagegroups[imagegroup_id] = imagegroup
    end
    return slide
end

function Experiment(experiment_path::String)
    parts = splitpath(experiment_path)
    exname = last(parts)
    base_path = length(parts) == 1 ? "./" : joinpath(parts[1:end-1])
    Experiment(exname, base_path)
end

function Experiment(exname::String, base_path::String)
    experiment_path = joinpath(base_path, exname)
    slnames = readdir(experiment_path)
    experiment = Experiment{Slide}(exname, base_path)
    @info slnames
    for slname in slnames
        slide_files = _slide_file_names(experiment_path, slname)
        if !all(isfile, values(slide_files))
            @warn "Some required files missing for $slname, skipping"
            continue
        end
        slide = Slide(slname, experiment)
        experiment[slname] = slide
    end
    return experiment
end


#############
# Accessors #
#############

"""
    experiment_name(obj)

Return the experiment name for any `XXX` data object

- [`Experiment`](@ref)
- [`Slide`](@ref)
- [`ImageGroup`](@ref)
- [`FOV`](@ref)

For slides, imagegroups, or fovs, returns the parent experiment name.
See also [`experiment_path`](@ref), [`slide_name`](@ref).
"""
experiment_name(ex::Experiment) = ex.name
experiment_name(sl::Slide) = experiment_name(sl.experiment)
experiment_name(gr::ImageGroup) = experiment_name(gr.slide)
experiment_name(fov::FOV) = experiment_name(fov.imagegroup)


"""
    experiment_path(obj)

Return the relative path to the experiment for any `XXX` data object

- [`Experiment`](@ref)
- [`Slide`](@ref)
- [`ImageGroup`](@ref)
- [`FOV`](@ref)

For slides, imagegroups, or fovs, returns the parent experiment path.
See also [`slide_path`](@ref), [`imagegroups_path`](@ref).
"""
experiment_path(ex::Experiment) = joinpath(ex.base_path, ex.name)
experiment_path(sl::Slide) = experiment_path(sl.experiment)
experiment_path(gr::ImageGroup) = experiment_path(gr.slide)
experiment_path(fov::FOV) = experiment_path(fov.imagegroup)

"""
    slide_path(obj)

Return the relative path to the slide(s) for any `XXX` data object,

- [`Experiment`](@ref)
- [`Slide`](@ref)
- [`ImageGroup`](@ref)
- [`FOV`](@ref)

For experiments, returns a vector of paths for all slides.
For imagegroups, or fovs, returns the parent slide path.
See also [`slide_path`](@ref), [`imagegroups_path`](@ref).
"""
slide_path(ex::Experiment) = [slide_path(slide) for slide in values(ex.slides)]
slide_path(sl::Slide) = joinpath(experiment_path(sl), sl.name)
slide_path(gr::ImageGroup) = slide_path(gr.slide)
slide_path(fov::FOV) = slide_path(fov.imagegroup)

"""
    imagegroups_path(obj)

Return the relative path to imagegroup images for any `XXX` data object

- [`Experiment`](@ref)
- [`Slide`](@ref)
- [`ImageGroup`](@ref)
- [`FOV`](@ref)

This is `{EXPERIMENT_ROOT}/{SLIDE_NAME}/imagegroups`.

For experiments, returns a vector of paths to imagegroups for all slides.
For imagegroups or fovs, returns the path to imagegroups for the parent slide.
See also [`fovs_path`](@ref), [`canonical_image_path`](@ref).
"""
imagegroups_path(ex::Experiment) = [imagegroups_path(slide) for slide in values(ex.slides)]
imagegroups_path(sl::Slide) = joinpath(slide_path(sl), "imagegroups")
imagegroups_path(gr::ImageGroup) = joinpath(slide_path(gr), "imagegroups")
imagegroups_path(fov::FOV) = imagegroup_path(fov.imagegroup)

"""
    fovs_path(obj)

Return the relative path to fov images for any `XXX` data object

- [`Experiment`](@ref)
- [`Slide`](@ref)
- [`ImageGroup`](@ref)
- [`FOV`](@ref)

This is `{EXPERIMENT_ROOT}/{SLIDE_NAME}/fovs`.

For experiments, returns a vector of paths to fovs for all slides.
For imagegroups or fovs, returns the path to fovs for the parent slide.
See also [`imagegroups_path`](@ref), [`canonical_image_path`](@ref).
"""
fovs_path(ex::Experiment) = [fovs_path(slide) for slide in values(ex.slides)]
fovs_path(sl::Slide) = joinpath(slide_path(sl), "fovs")
fovs_path(gr::ImageGroup) = joinpath(slide_path(gr), "fovs")
fovs_path(fov::FOV) = joinpath(slide_path(fov), "fovs")

"""
    slide_names(ex::Experiment)

Returns a vector of names for slides in `ex`.
Equivalent to `collect(keys(ex.slides))`
"""
slide_names(ex::Experiment) = collect(keys(ex.slides))

"""
    slide_name(obj)

Return the name of a [`Slide`],
or the parent slide of a [`ImageGroup`](@ref) or [`FOV`](@ref).
"""
slide_name(sl::Slide) = sl.name
slide_name(gr::ImageGroup) = slide_name(gr.slide)
slide_name(fov::FOV) = slide_name(fov.imagegroup)

"""
    imagegroup_name(obj)

Return the name of a [`ImageGroup`],
or the parent imagegroup of a [`FOV`](@ref).
"""
imagegroup_name(gr::ImageGroup) = gr.id
imagegroup_name(fov::FOV) = imagegroup_name(fov.imagegroup)

fov_ids(gr::ImageGroup) = [fov.id for fov in fovs(gr)]
fovs(gr::ImageGroup) = gr.fovs

"""
    transcripts(obj)

Return a transcripts table for any `XXX` data object
other than an experiment.

- [`Slide`](@ref)
- [`ImageGroup`](@ref)
- [`FOV`](@ref)

For slides and imagegroups, a `GroupedDataFrame` is returned,
grouped on fov indices,
meaning they can be indexed with a group key.

```julia-repl
julia> transcripts(imagegroup)[(; fov=24)]
141879×13 SubDataFrame
    Row │ fov     cell_ID  cell        x_local_px  y_local_px  x_global_px  y_global_px     z      target     CellComp   imagegroup_id  x_imagegroup_px  y_imagegroup_px
        │ UInt16  UInt16   String15    UInt16      UInt16      Float32      Float32         Int16  String31   String15   Int16    Int64      Int64
────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
      1 │     24        0  c_2_24_0          4254         499      64439.6       1.02615e5      8  Blk        None             5      12686       4715
      2 │     24        0  c_2_24_0          4256         738      64441.3       1.02376e5      7  Abl2       None             5      12688       4954
      3 │     24        0  c_2_24_0          4255         738      64441.0       1.02376e5      8  Abl2       None             5      12687       4954
   ⋮    │   ⋮        ⋮         ⋮           ⋮           ⋮            ⋮             ⋮           ⋮        ⋮          ⋮         ⋮         ⋮          ⋮
 141877 │     24      431  c_2_24_431         218        2728      60403.9       1.00386e5      7  C5ar2      Cytoplasm        5       8650       6944
 141878 │     24      431  c_2_24_431         212        2713      60397.4       1.00401e5      6  Cxcr6      Cytoplasm        5       8644       6929
 141879 │     24      431  c_2_24_431         210        2715      60395.8       1.00399e5      4  Fgf13      Cytoplasm        5       8642       6931
                                                                                                                                    141819 rows omitted
```
"""
transcripts(sl::Slide) = sl.transcripts
transcripts(gr::ImageGroup) = transcripts(gr.slide)[[(; fov) for fov in vec(fov_ids(gr)) if haskey(gr.slide, (; fov))]]
transcripts(fov::FOV) = transcripts(fov.imagegroup.slide)[(; fov=fov.id)]

positions(sl::Slide) = sl.positions
positions(gr::ImageGroup) = subset(positions(gr.slide), "imagegroup_id" => ByRow(==(imagegroup_name(gr))))
positions(fov::FOV) = subset(positions(fov.imagegroup.slide), "fov" => ByRow(==(fov.id)))

