
"""
    DiskStore(db)

Dict-like object that writes to disk
"""
struct DiskStore
    db::SQLite.DB
    function DiskStore(db::SQLite.DB)
        if isempty(DBInterface.execute(db, 
                "SELECT name FROM sqlite_master WHERE type='table' AND name='kv'")
                )
            DBInterface.execute(db, "CREATE TABLE kv(key TEXT, value BLOB)")
        end
        new(db)
    end
end

DiskStore(x::String) = DiskStore(SQLite.DB(x))

Base.keys(store::DiskStore) = (entry.key for entry in DBInterface.execute(store.db, "SELECT key FROM kv"))
Base.haskey(store::DiskStore, key) = key ∈ keys(store)

_query_on_key(store::DiskStore, key::String) = DBInterface.execute(store.db, "SELECT value from kv where key=?", (key,))

function Base.get(store::DiskStore, key::String, default)
    query = _query_on_key(store,key)
    if isempty(query)
        return default 
    else
        return first(query).value
    end
end

function Base.getindex(store::DiskStore, key::String)
    first(_query_on_key(store,value)).value
end

function Base.setindex!(store::DiskStore, value, key::String)
    DBInterface.execute(store.db, "DELETE FROM kv WHERE key=?", (key,))
    DBInterface.execute(store.db, "INSERT INTO kv VALUES (?, ?)", (key, value))
end

_properties_path(base_path, name) = joinpath(base_path, "$(name)_properties.sqlite")

function _build_or_get_props(base_path, name)
    propspath = _properties_path(base_path, name)
    !isdir(dirname(propspath)) && mkpath(dirname(propspath))
    props = DiskStore(propspath)
    haskey(props, "name") || setindex!(props, name, "name")
    return props
end
