# SpatialViz/src/backends/napari_backend.jl
# Napari backend via PythonCall.jl.

"""
    NapariBackend

Renders spatial data in a Napari viewer window via PythonCall.jl.
Requires a working Python environment with `napari` installed.
"""
struct NapariBackend <: VisualizationBackend end

function create_renderer(::Type{NapariBackend}, interactive::Bool)
    # TODO Phase 2: launch or reuse a Napari viewer via PythonCall
    # viewer = pyimport("napari").Viewer()
    error("NapariBackend.create_renderer: not yet implemented (Phase 2 deliverable)")
end
