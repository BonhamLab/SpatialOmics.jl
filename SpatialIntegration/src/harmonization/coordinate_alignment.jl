# SpatialIntegration/src/harmonization/coordinate_alignment.jl
# Spatial registration of multiple datasets to a common coordinate system.

"""
    align_coordinate_systems(datasets::Vector{SpatialDataset};
                              method::Symbol = :landmark,
                              reference_idx::Int = 1) -> Vector{SpatialDataset}

Register all datasets to the coordinate system of `datasets[reference_idx]`.
Returns a new vector of `SpatialDataset` with updated transformations; the
original datasets are unchanged.

Supported `method` values (planned):
- `:landmark`     — manual or detected landmark-based rigid registration
- `:icp`          — Iterative Closest Point alignment
- `:elastix`      — deformable registration via Elastix (requires PythonCall)
"""
function align_coordinate_systems(
    datasets::Vector{SpatialDataset};
    method::Symbol = :landmark,
    reference_idx::Int = 1,
)
    # TODO Phase 3: implement
    error("align_coordinate_systems: not yet implemented (Phase 3 deliverable)")
end
