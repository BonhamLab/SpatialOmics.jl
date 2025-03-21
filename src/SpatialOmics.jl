module SpatialOmics

export fov,
       fov!,
       fovs

using DataStructures
using DataFrames
using DimensionalData

using Images
using FileIO
using TOML
using Preferences

using GeometryOps
using GeometryBasics
import GeoInterface

import SQLite
import DBInterface

# Write your package code here.
include("datastores.jl")
include("data_structures.jl")
include("metadata.jl")
include("io.jl")

end
