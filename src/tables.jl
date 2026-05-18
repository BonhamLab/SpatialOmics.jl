# ── SpatialTable ──────────────────────────────────────────────────────────────

struct SpatialTable
    X            :: AbstractMatrix         # (n_obs × n_var)
    obs          :: Any                    # NamedTuple of Vectors — per-obs metadata
    var          :: Any                    # NamedTuple of Vectors — per-var metadata
    region       :: Union{String, Nothing} # element name this table annotates
    region_key   :: Symbol
    instance_key :: Symbol                 # which obs column links to element instance_id
    region_kind  :: Symbol                 # :shapes | :labels | :points | :none
end

function SpatialTable(X::AbstractMatrix;
                      obs          = NamedTuple(),
                      var          = NamedTuple(),
                      region       = nothing,
                      region_key   = :region,
                      instance_key = :instance_id,
                      region_kind  = :none)
    SpatialTable(X, obs, var, region, region_key, instance_key, region_kind)
end

nobs(tbl::SpatialTable)      = size(tbl.X, 1)
nvar(tbl::SpatialTable)      = size(tbl.X, 2)

function var_names(tbl::SpatialTable)
    tbl.var isa NamedTuple && haskey(tbl.var, :name) ?
        tbl.var.name : string.(1:nvar(tbl))
end

function _obs_instance_ids(tbl::SpatialTable)
    key = tbl.instance_key
    tbl.obs isa NamedTuple && haskey(tbl.obs, key) ?
        tbl.obs[key] : Int32.(1:nobs(tbl))
end

# ── Dataset typed accessor ─────────────────────────────────────────────────────

function tables(ds::SpatialDataset, name::String)
    el = ds.elements[name]
    el isa SpatialTable || error("Element \"$name\" is not SpatialTable (got $(typeof(el)))")
    el
end

# ── feature API — full dataset (all obs in shape order) ───────────────────────

function feature(ds::SpatialDataset, gene::String; region::String)
    shape_el = ds.elements[region]
    shape_el isa SpatialShapes || error("Element \"$region\" is not SpatialShapes")
    ids_ordered = shape_el.instance_id

    for (_, el) in ds.elements
        el isa SpatialTable || continue
        el.region == region || continue
        idx = findfirst(==(gene), var_names(el))
        idx === nothing && error("Gene \"$gene\" not found in table linked to \"$region\"")
        obs_ids   = _obs_instance_ids(el)
        id_to_row = Dict(id => i for (i, id) in enumerate(obs_ids))
        return [el.X[id_to_row[id], idx] for id in ids_ordered]
    end
    error("No SpatialTable linked to element \"$region\"")
end
