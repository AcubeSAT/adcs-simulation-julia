struct Optional{T}
    data::Union{Nothing,T}
    Optional{T}(data::Union{Nothing,T} = nothing) where {T} = new{T}(data)
end

Base.getindex(s::Optional) = s.data
function Base.isequal(a::Optional, b::Optional)
    if a.data === nothing || b.data === nothing
        return true
    else
        return isequal(a.data, b.data)
    end
end

Base.convert(::Type{Optional{T}}, ::Nothing) where {T} = Optional{T}()
Base.convert(::Type{Optional{T}}, x) where {T} = Optional{T}(convert(T, x))
Base.convert(::Type{Optional{T}}, x::Optional) where {T} = convert(Optional{T}, x[])
