"""
    CellInteraction

Cell-cell communication and neighbourhood analysis for spatial datasets.
Provides ligand-receptor analysis, spatial neighbourhood profiling, and
communication scoring.
"""
module CellInteraction

using Statistics
using ThreadsX
using SpatialOmicsBase

export ligand_receptor_analysis, neighborhood_analysis, communication_score

"""
    ligand_receptor_analysis(ds::SpatialDataset;
                             database::Symbol = :cellchat,
                             method::Symbol = :proximity_ligation,
                             radius::Real = 50.0) -> DataFrame

Screen for active ligand-receptor pairs between spatially-adjacent cell types.
Returns a DataFrame with columns: `ligand`, `receptor`, `sender_type`,
`receiver_type`, `score`, `pvalue`.

Supported `database` values (planned): `:cellchat`, `:nichenet`, `:liana`
Supported `method` values (planned):
- `:proximity_ligation` — score pairs based on spatial co-localisation
- `:expression_only`    — expression-based scoring without spatial constraint
"""
function ligand_receptor_analysis(
    ds::SpatialDataset;
    database::Symbol = :cellchat,
    method::Symbol = :proximity_ligation,
    radius::Real = 50.0,
)
    # TODO Phase 3: implement
    error("ligand_receptor_analysis: not yet implemented (Phase 3 deliverable)")
end

"""
    neighborhood_analysis(ds::SpatialDataset;
                          cell_type_key::String = "cell_type",
                          radius::Real = 100.0,
                          n_permutations::Int = 1000) -> DataFrame

Profile the spatial composition of cell-type neighbourhoods. Returns a
cell-type × cell-type DataFrame of enrichment scores (observed / expected).
"""
function neighborhood_analysis(
    ds::SpatialDataset;
    cell_type_key::String = "cell_type",
    radius::Real = 100.0,
    n_permutations::Int = 1000,
)
    # TODO Phase 3: implement; parallelise permutations with ThreadsX
    error("neighborhood_analysis: not yet implemented (Phase 3 deliverable)")
end

"""
    communication_score(ds::SpatialDataset,
                        sender_type::String,
                        receiver_type::String;
                        method::Symbol = :cellchat) -> Float64

Compute an aggregate cell-cell communication score between `sender_type` and
`receiver_type` populations in `ds`.
"""
function communication_score(
    ds::SpatialDataset,
    sender_type::String,
    receiver_type::String;
    method::Symbol = :cellchat,
)
    # TODO Phase 3: implement
    error("communication_score: not yet implemented (Phase 3 deliverable)")
end

end # module CellInteraction
