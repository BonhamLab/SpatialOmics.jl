# ── Pretty REPL display ────────────────────────────────────────────────────────

function Base.show(io::IO, fmt::CosMx)
    m = fmt.morphology_dir === nothing ? "" : "morphology_dir=$(repr(fmt.morphology_dir))"
    print(io, "CosMx($m)")
end

function Base.show(io::IO, cs::CoordinateSystem)
    print(io, "CoordinateSystem(\"$(cs.name)\", $(cs.axes[1])/$(cs.axes[2]), $(cs.units[1])/$(cs.units[2]))")
end

Base.show(io::IO, t::Identity)  = print(io, "Identity: \"$(t.src)\" → \"$(t.dst)\"")
Base.show(io::IO, t::Affine)    = print(io, "Affine: \"$(t.src)\" → \"$(t.dst)\"")
Base.show(io::IO, t::Sequence)  = print(io, "Sequence($(length(t.steps)) steps): \"$(t.src)\" → \"$(t.dst)\"")

function Base.show(io::IO, pts::SpatialPoints{T}) where T
    n  = length(pts)
    nf = length(pts.feature_codebook)
    fstr = nf == 0   ? "no features" :
           nf <= 5   ? join(pts.feature_codebook, ", ") :
                       join(pts.feature_codebook[1:5], ", ") * " … ($nf total)"
    cs = isempty(pts.coord_system) ? "" : " [$(pts.coord_system)]"
    print(io, "SpatialPoints{$T}: $n points · $fstr$cs")
end

function Base.show(io::IO, shp::SpatialShapes{G}) where G
    cs = isempty(shp.coord_system) ? "" : " [$(shp.coord_system)]"
    print(io, "SpatialShapes{$(nameof(G))}: $(length(shp)) shapes$cs")
end

function Base.show(io::IO, ext::SpatialExtent)
    cs = isempty(ext.coord_system) ? "" : " [$(ext.coord_system)]"
    print(io, "SpatialExtent(x=$(ext.xmin)..$(ext.xmax), y=$(ext.ymin)..$(ext.ymax))$cs")
end

function Base.show(io::IO, roi::SpatialROI{G}) where G
    ext = roi.extent
    cs  = isempty(roi.coord_system) ? "" : " [$(roi.coord_system)]"
    print(io, "SpatialROI{$(nameof(G))}(x=$(ext.xmin)..$(ext.xmax), y=$(ext.ymin)..$(ext.ymax))$cs")
end

function Base.show(io::IO, lbl::SpatialLabels{T}) where T
    sz = join(size(lbl.data), "×")
    ni = length(lbl.instance_map)
    cs = isempty(lbl.coord_system) ? "" : " [$(lbl.coord_system)]"
    print(io, "SpatialLabels{$T}: $sz, $ni instances$cs")
end

function Base.show(io::IO, tbl::SpatialTable)
    reg = tbl.region === nothing ? "unlinked" : "→ \"$(tbl.region)\""
    print(io, "SpatialTable: $(nobs(tbl)) obs × $(nvar(tbl)) var ($reg)")
end

function Base.show(io::IO, img::SpatialImage{T}) where T
    sz  = join(size(img.data), "×")
    nc  = nchannels(img)
    nl  = length(img.pyramid)
    cs  = isempty(img.coord_system) ? "" : " [$(img.coord_system)]"
    pyr = nl > 0 ? ", $nl pyramid level$(nl == 1 ? "" : "s")" : ""
    print(io, "SpatialImage{$T}: $sz$pyr$cs")
end

function Base.show(io::IO, v::SpatialElementView)
    mode = v.overlap == :full ? " overlap=:full" : ""
    print(io, "SpatialElementView{$(typeof(v.parent))}$mode (lazy) ↩ ")
    show(io, v.roi)
end

function Base.show(io::IO, v::SpatialDatasetView)
    n = length(v.parent.elements)
    print(io, "SpatialDatasetView($n element$(n == 1 ? "" : "s"), lazy) ↩ ")
    show(io, v.roi)
end

# Compact inline form — used when ds appears as a field or inside a larger structure
function Base.show(io::IO, ds::SpatialDataset)
    n   = length(ds.elements)
    ncs = length(ds.coord_systems)
    cs  = ncs == 0 ? "" :
          ncs <= 3  ? " [$(join(keys(ds.coord_systems), ", "))]" :
                      " ($ncs coord systems)"
    print(io, "SpatialDataset($n element$(n == 1 ? "" : "s")$cs)")
end

# Full REPL form — used when ds is displayed at top level
function Base.show(io::IO, ::MIME"text/plain", ds::SpatialDataset)
    n   = length(ds.elements)
    ncs = length(ds.coord_systems)
    println(io, "SpatialDataset with $n element$(n == 1 ? "" : "s") and $ncs coord_system$(ncs == 1 ? "" : "s"):")
    for (name, el) in ds.elements
        print(io, "  \"$name\" => ")
        show(io, el)
        println(io)
    end
    if !isempty(ds.coord_systems)
        css  = collect(values(ds.coord_systems))
        ncs  = length(css)
        head = join(["\"$(cs.name)\" ($(cs.units[1]))" for cs in css[1:min(3,ncs)]], ", ")
        tail = ncs > 3 ? " … ($ncs total)" : ""
        print(io, "  coord_systems: ", head, tail)
        if !isempty(ds.transforms)
            print(io, "\n  transforms: $(length(ds.transforms))")
        end
    end
end
