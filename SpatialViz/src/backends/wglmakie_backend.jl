# SpatialViz/src/backends/wglmakie_backend.jl
# WGLMakie + Bonito backend for browser-based / Jupyter deployment.

"""
    WGLMakieBackend

Renders spatial data in the browser using WGLMakie.jl. Can be served via
Bonito.jl for standalone web deployment or used inside Jupyter/Pluto notebooks.
"""
struct WGLMakieBackend <: VisualizationBackend end

"""
    BonitoBackend

Alias for `WGLMakieBackend` that emphasises Bonito.jl web-server deployment.
"""
const BonitoBackend = WGLMakieBackend

function create_renderer(::Type{WGLMakieBackend}, interactive::Bool)
    # TODO Phase 2: return a WGLMakie Figure (activates WGLMakie backend)
    error("WGLMakieBackend.create_renderer: not yet implemented (Phase 2 deliverable)")
end
