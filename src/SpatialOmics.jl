module SpatialOmics

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
include("metadata.jl")
include("data_structures.jl")
include("io.jl")

end
