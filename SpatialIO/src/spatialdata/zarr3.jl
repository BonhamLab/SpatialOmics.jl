# SpatialIO/src/spatialdata/zarr3.jl
#
# Low-level utilities for reading and writing Zarr v3 stores as produced by
# spatialdata (Python).  The format differs from Zarr v2 in several ways:
#   - Metadata is in `zarr.json` (not `.zarray`/`.zattrs`)
#   - Chunk paths are prefixed with `c/`: chunk (i,j,k) lives at `c/i/j/k`
#   - Supported codecs: `bytes` (little-endian raw) + `zstd`, and
#     `vlen-utf8` (variable-length strings) + `zstd`
#
# This module is internal to SpatialIO.

using CodecZstd: ZstdDecompressor, ZstdCompressor
using JSON3
using DiskArrays

# ─── Metadata helpers ─────────────────────────────────────────────────────────

"""Read and parse `zarr.json` from `path`."""
_zread_meta(path::String) = JSON3.read(Base.read(joinpath(path, "zarr.json"), String))

"""
List immediate child groups (sub-directories that contain a `zarr.json`)
under `path`, in sorted order.
"""
function _zlist_groups(path::String)::Vector{String}
    isdir(path) || return String[]
    sort([e for e in readdir(path)
          if isdir(joinpath(path, e)) && isfile(joinpath(path, e, "zarr.json"))])
end

# ─── dtype ────────────────────────────────────────────────────────────────────

function _zdtype_to_julia(dt::String)
    dt == "bool"    && return Bool
    dt == "uint8"   && return UInt8
    dt == "uint16"  && return UInt16
    dt == "uint32"  && return UInt32
    dt == "uint64"  && return UInt64
    dt == "int8"    && return Int8
    dt == "int16"   && return Int16
    dt == "int32"   && return Int32
    dt == "int64"   && return Int64
    dt == "float32" && return Float32
    dt == "float64" && return Float64
    error("_zdtype_to_julia: unknown dtype \"$dt\"")
end

function _zjulia_to_dtype(::Type{T}) where T
    T == Bool    && return "bool"
    T == UInt8   && return "uint8"
    T == UInt16  && return "uint16"
    T == UInt32  && return "uint32"
    T == UInt64  && return "uint64"
    T == Int8    && return "int8"
    T == Int16   && return "int16"
    T == Int32   && return "int32"
    T == Int64   && return "int64"
    T == Float32 && return "float32"
    T == Float64 && return "float64"
    error("_zjulia_to_dtype: unsupported type $T")
end

# ─── Compression ──────────────────────────────────────────────────────────────

function _zdecompress(raw::Vector{UInt8}, codecs)::Vector{UInt8}
    for c in codecs
        String(c.name) == "zstd" && return transcode(ZstdDecompressor, raw)
    end
    return raw
end

function _zcompress(data::Vector{UInt8}; level::Int = 0)::Vector{UInt8}
    transcode(ZstdCompressor(; level = level), data)
end

# ─── vlen-utf8 encode / decode ────────────────────────────────────────────────
#
# Format (numcodecs VLenUTF8, used by zarr-python):
#   [N: uint32 LE] followed by N × ([len_i: uint32 LE] [utf8 bytes…])

function _decode_vlen_utf8(bytes::Vector{UInt8})::Vector{String}
    n = Int(reinterpret(UInt32, @view bytes[1:4])[1])
    result = Vector{String}(undef, n)
    pos = 5
    for i in 1:n
        len = Int(reinterpret(UInt32, @view bytes[pos:pos+3])[1])
        pos += 4
        result[i] = String(copy(bytes[pos:pos+len-1]))
        pos += len
    end
    return result
end

function _encode_vlen_utf8(strs::Vector{String})::Vector{UInt8}
    n       = length(strs)
    encoded = [Vector{UInt8}(s) for s in strs]
    total   = 4 + sum(4 + length(b) for b in encoded; init=0)
    buf     = Vector{UInt8}(undef, total)
    copyto!(buf, 1, reinterpret(UInt8, [UInt32(n)]), 1, 4)
    pos = 5
    for b in encoded
        copyto!(buf, pos, reinterpret(UInt8, [UInt32(length(b))]), 1, 4)
        pos += 4
        copyto!(buf, pos, b, 1, length(b))
        pos += length(b)
    end
    return buf
end

# ─── Chunk path ───────────────────────────────────────────────────────────────
#
# Zarr v3 "default" chunk key encoding with "/" separator:
# chunk at 0-based grid coordinates (c0, c1, …, cN) → "c/c0/c1/…/cN"

function _chunk_path(base::String, coords)
    joinpath(base, "c", join(string.(coords), "/"))
end

# ─── Eager 1-D array reader ──────────────────────────────────────────────────
#
# Used for obs/var column arrays (always 1-D) and for sparse-matrix components.

function _zread_array_1d(path::String)
    meta   = _zread_meta(path)
    shape  = Int.(meta.shape)
    length(shape) == 1 || error("_zread_array_1d: expected 1-D, got shape $shape at $path")
    codecs = meta.codecs
    cs     = Int(meta.chunk_grid.configuration.chunk_shape[1])
    n      = shape[1]
    n_ch   = cld(n, cs)

    if any(c -> String(c.name) == "vlen-utf8", codecs)
        result = String[]
        sizehint!(result, n)
        for i in 0:n_ch-1
            cpath = _chunk_path(path, (i,))
            isfile(cpath) || continue
            chunk_strs = _decode_vlen_utf8(_zdecompress(Base.read(cpath), codecs))
            # Last chunk may be padded to full chunk size; truncate to actual remaining
            remaining = n - length(result)
            append!(result, view(chunk_strs, 1:min(length(chunk_strs), remaining)))
        end
        return result
    else
        T      = _zdtype_to_julia(String(meta.data_type))
        result = Vector{T}(undef, n)
        for i in 0:n_ch-1
            c_start = i * cs + 1
            c_end   = min(c_start + cs - 1, n)
            cpath   = _chunk_path(path, (i,))
            if isfile(cpath)
                dec = _zdecompress(Base.read(cpath), codecs)
                vals = reinterpret(T, dec)
                result[c_start:c_end] .= vals[1:(c_end - c_start + 1)]
            else
                result[c_start:c_end] .= zero(T)
            end
        end
        return result
    end
end

# ─── Eager N-D array reader (small arrays) ───────────────────────────────────

function _zread_array_nd(path::String)
    meta        = _zread_meta(path)
    shape       = Int.(meta.shape)
    N           = length(shape)
    chunk_shape = Int.(meta.chunk_grid.configuration.chunk_shape)
    codecs      = meta.codecs
    T           = _zdtype_to_julia(String(meta.data_type))

    result = Array{T}(undef, shape...)
    n_chunks = ntuple(i -> cld(shape[i], chunk_shape[i]), N)

    for chunk_coords in Iterators.product(ntuple(i -> 0:n_chunks[i]-1, N)...)
        g_start = ntuple(i -> chunk_coords[i] * chunk_shape[i] + 1, N)
        g_end   = ntuple(i -> min((chunk_coords[i]+1) * chunk_shape[i], shape[i]), N)
        c_sizes = ntuple(i -> g_end[i] - g_start[i] + 1, N)
        out_r   = ntuple(i -> g_start[i]:g_end[i], N)

        cpath = _chunk_path(path, chunk_coords)
        if isfile(cpath)
            dec        = _zdecompress(Base.read(cpath), codecs)
            chunk_c    = reshape(reinterpret(T, dec), reverse(chunk_shape))
            chunk_full = permutedims(chunk_c, N:-1:1)
            local_r    = ntuple(i -> 1:c_sizes[i], N)
            result[out_r...] .= chunk_full[local_r...]
        else
            result[out_r...] .= zero(T)
        end
    end
    return result
end

# ─── ZarrV3Array — lazy DiskArray for images and labels ──────────────────────

struct ZarrV3Array{T,N} <: DiskArrays.AbstractDiskArray{T,N}
    path        ::String
    shape       ::NTuple{N,Int}
    chunk_shape ::NTuple{N,Int}
    codecs      ::Any           # JSON3 array of codec objects
    fill_value  ::T
end

function ZarrV3Array(path::String)
    meta = _zread_meta(path)
    String(meta.node_type) == "array" ||
        error("ZarrV3Array: not an array node at $path")
    T  = _zdtype_to_julia(String(meta.data_type))
    N  = length(meta.shape)
    shape       = NTuple{N,Int}(Int.(meta.shape))
    chunk_shape = NTuple{N,Int}(Int.(meta.chunk_grid.configuration.chunk_shape))
    fv_raw = meta.fill_value
    fv     = isa(fv_raw, Number) ? T(fv_raw) : zero(T)
    return ZarrV3Array{T,N}(path, shape, chunk_shape, meta.codecs, fv)
end

Base.size(a::ZarrV3Array) = a.shape
DiskArrays.haschunks(::ZarrV3Array) = DiskArrays.Chunked()
DiskArrays.eachchunk(a::ZarrV3Array) = DiskArrays.GridChunks(a, a.chunk_shape)

function DiskArrays.readblock!(
    a    ::ZarrV3Array{T,N},
    aout ::AbstractArray,
    r    ::OrdinalRange...,
) where {T,N}
    cs    = a.chunk_shape
    shape = a.shape

    ci_start = ntuple(i -> div(first(r[i]) - 1, cs[i]), N)
    ci_end   = ntuple(i -> div(last(r[i])  - 1, cs[i]), N)

    for chunk_coords in Iterators.product(ntuple(i -> ci_start[i]:ci_end[i], N)...)
        g_start = ntuple(i -> chunk_coords[i] * cs[i] + 1, N)
        g_end   = ntuple(i -> min((chunk_coords[i]+1) * cs[i], shape[i]), N)

        # Intersection of this chunk with the requested range
        inter_start = ntuple(i -> max(Int(first(r[i])), g_start[i]), N)
        inter_end   = ntuple(i -> min(Int(last(r[i])),  g_end[i]),   N)
        any(i -> inter_start[i] > inter_end[i], 1:N) && continue

        local_r = ntuple(i -> (inter_start[i] - g_start[i] + 1):(inter_end[i] - g_start[i] + 1), N)
        aout_r  = ntuple(i -> (inter_start[i] - Int(first(r[i])) + 1):(inter_end[i] - Int(first(r[i])) + 1), N)

        cpath = _chunk_path(a.path, chunk_coords)
        if isfile(cpath)
            dec = _zdecompress(Base.read(cpath), a.codecs)
            # Zarr v3 stores chunks in C order (row-major, last axis varies fastest).
            # Julia uses Fortran order (column-major, first axis varies fastest).
            # Reshape with reversed dims then permute back to get correct indexing.
            chunk_c    = reshape(reinterpret(T, dec), reverse(cs))
            chunk_data = permutedims(chunk_c, N:-1:1)
            aout[aout_r...] .= chunk_data[local_r...]
        else
            aout[aout_r...] .= a.fill_value
        end
    end
    return aout
end

# ─── Array writer ─────────────────────────────────────────────────────────────

"""Write a Julia `AbstractArray` as a Zarr v3 array at `path` with zstd compression."""
function _zwrite_array(
    path        ::String,
    data        ::AbstractArray{T,N};
    chunk_shape ::NTuple{N,Int} = size(data),  # default: one chunk
    zstd_level  ::Int = 0,
) where {T,N}
    mkpath(joinpath(path, "c"))

    shape     = size(data)
    dtype_str = _zjulia_to_dtype(T)
    n_chunks  = ntuple(i -> cld(shape[i], chunk_shape[i]), N)

    # Write zarr.json
    meta = Dict(
        "zarr_format" => 3,
        "node_type"   => "array",
        "shape"       => collect(shape),
        "data_type"   => dtype_str,
        "chunk_grid"  => Dict(
            "name"          => "regular",
            "configuration" => Dict("chunk_shape" => collect(chunk_shape)),
        ),
        "chunk_key_encoding" => Dict(
            "name"          => "default",
            "configuration" => Dict("separator" => "/"),
        ),
        "fill_value" => 0,
        "codecs"     => [
            Dict("name" => "bytes", "configuration" => Dict("endian" => "little")),
            Dict("name" => "zstd",  "configuration" => Dict("level" => zstd_level, "checksum" => false)),
        ],
        "attributes"          => Dict{String,Any}(),
        "storage_transformers" => [],
    )
    Base.write(joinpath(path, "zarr.json"), JSON3.write(meta))

    # Write chunks
    for chunk_coords in Iterators.product(ntuple(i -> 0:n_chunks[i]-1, N)...)
        g_start = ntuple(i -> chunk_coords[i] * chunk_shape[i] + 1, N)
        g_end   = ntuple(i -> min((chunk_coords[i]+1) * chunk_shape[i], shape[i]), N)
        slice   = ntuple(i -> g_start[i]:g_end[i], N)
        chunk   = collect(data[slice...])
        # Write in zarr C order (last index varies fastest): permute dims then vec in Fortran order.
        raw     = reinterpret(UInt8, vec(permutedims(chunk, N:-1:1)))
        compressed = _zcompress(collect(raw); level = zstd_level)
        cpath = _chunk_path(path, chunk_coords)
        mkpath(dirname(cpath))
        Base.write(cpath, compressed)
    end
end

"""Write a `Vector{String}` as a vlen-utf8 Zarr v3 array."""
function _zwrite_string_array(
    path       ::String,
    strs       ::Vector{String};
    zstd_level ::Int = 0,
)
    n = length(strs)
    mkpath(joinpath(path, "c"))

    meta = Dict(
        "zarr_format" => 3,
        "node_type"   => "array",
        "shape"       => [n],
        "data_type"   => "string",
        "chunk_grid"  => Dict(
            "name"          => "regular",
            "configuration" => Dict("chunk_shape" => [n]),  # one chunk
        ),
        "chunk_key_encoding" => Dict(
            "name"          => "default",
            "configuration" => Dict("separator" => "/"),
        ),
        "fill_value" => "",
        "codecs"     => [
            Dict("name" => "vlen-utf8", "configuration" => Dict{String,Any}()),
            Dict("name" => "zstd",      "configuration" => Dict("level" => zstd_level, "checksum" => false)),
        ],
        "attributes"          => Dict{String,Any}(),
        "storage_transformers" => [],
    )
    Base.write(joinpath(path, "zarr.json"), JSON3.write(meta))

    raw        = _encode_vlen_utf8(strs)
    compressed = _zcompress(raw; level = zstd_level)
    cpath      = _chunk_path(path, (0,))
    mkpath(dirname(cpath))
    Base.write(cpath, compressed)
end

"""Write a group zarr.json with the given attributes (Dict or JSON3.Object)."""
function _zwrite_group(path::String, attrs = Dict{String,Any}())
    mkpath(path)
    # Convert JSON3.Object → Dict so JSON3.write round-trips cleanly
    attrs_dict = _to_plain_dict(attrs)
    meta = Dict(
        "zarr_format" => 3,
        "node_type"   => "group",
        "attributes"  => attrs_dict,
    )
    Base.write(joinpath(path, "zarr.json"), JSON3.write(meta))
end

# Recursively convert JSON3.Object/Array to plain Julia dicts/arrays
function _to_plain_dict(x::JSON3.Object)
    Dict{String,Any}(String(k) => _to_plain_dict(v) for (k, v) in x)
end
function _to_plain_dict(x::JSON3.Array)
    [_to_plain_dict(v) for v in x]
end
_to_plain_dict(x) = x
