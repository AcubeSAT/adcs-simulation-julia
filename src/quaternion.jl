# Inspired by Michael Boyle and Base.Complex

# Signify standard arithmetic operations should be implemented on it
struct Quaternion{T<:Real} <: Number
    coeffs::SVector{4,T}
    Quaternion{T}(v::SVector{4,T}) where {T<:Real} = new{T}(v)
    # Always use SVector internally
    Quaternion{T}(v::A) where {T<:Real, A<:AbstractVector} = new{T}(SVector{4,T}(v))
end

# Construct without having to specify T
Quaternion(v::SVector{4,T}) where {T<:Real} = Quaternion{T}(v)
Quaternion(v::AbstractVector{T}) where {T<:Real} = Quaternion{T}(v)

# Explicitly state type, use SVector promote_rules for mixed types
(::Type{Quaternion{T}})(w, x, y, z) where {T<:Real} = Quaternion(SVector{4,T}(w, x, y, z))
# Rely on type inference, use SVector promote_rules for mixed types
function (::Type{Quaternion})(w, x, y, z)
    v = SVector{4}(w, x, y, z)
    return Quaternion{eltype(v)}(v)
end

coeffs(q::Quaternion) = getfield(q, :coeffs)

# Type-preserving copy constructor
Quaternion{T}(q::Quaternion{T}) where {T<:Real} = Quaternion(coeffs(q)...)
# Type conversion copy constructor
Quaternion{T}(q::Quaternion{S}) where {T<:Real, S<:Real} = Quaternion(SVector{4, T}(coeffs(q)...))

const QuaternionF64 = Quaternion{Float64}
const QuaternionF32 = Quaternion{Float32}
const QuaternionF16 = Quaternion{Float16}

# Given a numeric type, return a Quaternion (not instance) specialized on T
(::Type{Quaternion})(::Type{T}) where {T<:Real} = Quaternion{T}
# Given a Quaternion specialized on T, return another Quaternion also specialized on T
(::Type{Quaternion})(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}
