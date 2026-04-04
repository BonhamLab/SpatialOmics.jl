# SpatialViz/src/interactive/annotation.jl
# Manual cell-type annotation tools.

"""
    annotation_tool(ds::SpatialDataset;
                    cell_types::Vector{String} = String[],
                    backend::Type{<:VisualizationBackend} = MakieBackend,
                    save_path::Union{Nothing,String} = nothing)

Launch an interactive annotation session for `ds`. Users can paint regions
with cell-type labels, undo/redo, and save annotations to `save_path`.
Returns the annotated `SpatialDataset` when the session is closed.
"""
function annotation_tool(
    ds::SpatialDataset;
    cell_types::Vector{String} = String[],
    backend::Type{<:VisualizationBackend} = MakieBackend,
    save_path::Union{Nothing,String} = nothing,
)
    # TODO Phase 2: implement interactive annotation
    error("annotation_tool: not yet implemented (Phase 2 deliverable)")
end

"""
    save_annotations(ds::SpatialDataset, path::String)

Persist annotation labels from `ds.tables["annotations"]` to a CSV or JSON
file at `path`.
"""
function save_annotations(ds::SpatialDataset, path::String)
    error("save_annotations: not yet implemented (Phase 2 deliverable)")
end

"""
    load_annotations(ds::SpatialDataset, path::String) -> SpatialDataset

Load previously saved annotations from `path` into `ds` and return the
updated dataset.
"""
function load_annotations(ds::SpatialDataset, path::String)
    error("load_annotations: not yet implemented (Phase 2 deliverable)")
end
