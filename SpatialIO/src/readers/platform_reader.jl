# SpatialIO/src/readers/platform_reader.jl
# Abstract PlatformReader type and the universal load_spatial_data dispatcher.

"""
    PlatformReader

Abstract supertype for all vendor-format readers. Concrete subtypes implement
`read_data(reader, path)` and optionally `validate_path(reader, path)`.

New platforms are added by:
  1. Defining `struct MyPlatformReader <: PlatformReader end`
  2. Implementing `read_data(::MyPlatformReader, path::String)::SpatialDataset`
"""
abstract type PlatformReader end

# ---------------------------------------------------------------------------
# Canonical column schema for transcript SpatialPoints
#
# All platform readers MUST normalize to these names before constructing
# SpatialPoints so that downstream analysis code is platform-agnostic.
#
#   gene             — gene / target name (String)
#   cell_id          — integer cell instance identifier (Int32); links to
#                      SpatialShapes and SpatialTable via instance_key
#   cell             — string cell label, if the platform provides one
#   fov              — field-of-view index (Int32), if applicable
#   z                — z-plane index or depth (Int8 / Float32), if applicable
#   cell_compartment — subcellular compartment label (e.g. "Nuclear",
#                      "Cytoplasm"), if the platform provides one
#   quality_value    — per-transcript decode quality / confidence score
#
# Coordinate columns (x, y) are extracted into the coordinates matrix and
# do NOT appear in features.
# ---------------------------------------------------------------------------

"""
    read_data(reader::PlatformReader, path::String) -> SpatialDataset

Load a `SpatialDataset` from `path` using the given reader.
Must be implemented by every concrete `PlatformReader` subtype.
"""
function read_data end

"""
    validate_path(reader::PlatformReader, path::String) -> Bool

Check that `path` contains the expected files/structure for `reader`.
Default implementation returns `true`; override for stricter validation.
"""
validate_path(::PlatformReader, ::String) = true

"""
    load(reader::PlatformReader, path::String; validate::Bool=true) -> SpatialDataset

Load a `SpatialDataset` from `path` using the given reader instance. Reader
configuration is set at construction time via keyword arguments:

```julia
import SpatialIO as SIO
ds = SIO.load(SIO.CosMxReader(),                      "/data/my_experiment")
ds = SIO.load(SIO.CosMxReader(; morphology_dir=dir),  "/data/my_experiment")
```
"""
function load(reader::PlatformReader, path::String; validate::Bool = true)
    if validate && !validate_path(reader, path)
        error("Path does not look like valid $(typeof(reader)) output: $path")
    end
    return read_data(reader, path)
end

"""
    load_spatial_data(path::String, platform::Type{<:PlatformReader}; kwargs...)
        -> SpatialDataset

Internal dispatcher: construct a reader of type `platform`, optionally
validate `path`, then call `read_data`. Prefer `load(reader, path)` for new code.
"""
function load_spatial_data(
    path::String,
    platform::Type{<:PlatformReader};
    validate::Bool = true,
    kwargs...,
)
    reader = platform(; kwargs...)
    if validate && !validate_path(reader, path)
        error("Path does not look like a valid $(platform) output: $path")
    end
    return read_data(reader, path)
end
