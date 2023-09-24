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
Quaternion{T}(w, x, y, z) where {T<:Real} = Quaternion(SVector{4,T}(w, x, y, z))
# Rely on type inference, use SVector promote_rules for mixed types
function Quaternion(w, x, y, z)
    v = SVector{4}(w, x, y, z)
    return Quaternion{eltype(v)}(v)
end

# Type-preserving copy constructor
Quaternion{T}(Q::Quaternion{T}) where {T<:Real} = Quaternion(Q.coeffs...)
# Type conversion copy constructor
Quaternion{T}(Q::Quaternion{S}) where {T<:Real, S<:Real} = Quaternion(SVector{4, T}(Q.coeffs...))

const QuaternionF64 = Quaternion{Float64}
const QuaternionF32 = Quaternion{Float32}
const QuaternionF16 = Quaternion{Float16}

Base.zero(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}(zero(T), zero(T), zero(T), zero(T))
Base.zero(Q::Quaternion{T}) where {T<:Real} = Base.zero(typeof(Q))

Base.one(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}(one(T), zero(T), zero(T), zero(T))
Base.one(Q::Quaternion{T}) where {T<:Real} = Base.one(typeof(Q))

# Given a numeric type, return a Quaternion (not instance) specialized on T
(::Type{Quaternion})(::Type{T}) where {T<:Real} = Quaternion{T}
# Given a Quaternion specialized on T, return another Quaternion also specialized on T
(::Type{Quaternion})(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}

Base.eltype(::Type{Quaternion{T}}) where {T<:Real} = T
Base.promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T<:Real, S<:Real} = Quaternion{promote_type(T, S)}
Base.promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T<:Real, S<:Real} = Quaternion{promote_type(T, S)}

@inline function Base.getindex(Q::Quaternion, i::Integer)
    @boundscheck checkbounds(Q.coeffs, i)
    @inbounds return Q.coeffs[i]
end

# Let SVector handle the underlying indexing
# TODO: Should this return a view?
Base.@propagate_inbounds Base.getindex(Q::Quaternion, I) = @view Q.coeffs[I]

Base.length(::Quaternion) = 4

Base.real(::Type{Quaternion}) = eltype(Quaternion)
Base.real(::Type{Quaternion{T}}) where {T} = T
Base.real(Q::Quaternion{T}) where {T<:Real} = Q[1]
Base.imag(Q::Quaternion{T}) where {T<:Real} = @view Q.coeffs[2:4]
Base.vec(Q::Quaternion{T}) where {T<:Real} = @view Q.coeffs[2:4]

Base.widen(::Type{Quaternion{T}}) where {T} = Quaternion{widen(T)}
Base.big(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{big(T)}
Base.big(Q::Quaternion{T}) where {T<:Real} = Quaternion{big(T)}(Q)

Base.conj(Q::Quaternion) = Quaternion(Q[1], -Q[2], -Q[3], -Q[4])
Base.abs2(Q::Quaternion) = sum(abs2, Q.coeffs)
# TODO: Should do isnan/isinf checks and scale with max of abs of coeffs?
Base.abs(Q::Quaternion) = sqrt(abs2(Q))
Base.inv(Q::Quaternion) = conj(Q) / abs2(Q)
LinearAlgebra.norm(Q::Quaternion) = abs(Q)
LinearAlgebra.normalize(Q::Quaternion) = Q / abs(Q)

Base.:-(Q::Quaternion) = Quaternion(-Q.coeffs)

Base.:+(Q1::Quaternion, Q2::Quaternion) = Quaternion(Q1.coeffs + Q2.coeffs)
Base.:-(Q1::Quaternion, Q2::Quaternion) = Quaternion(Q1.coeffs - Q2.coeffs)
function Base.:*(Q1::Quaternion, Q2::Quaternion)
    return Quaternion(
        Q1[1]*Q2[1] - Q1[2]*Q2[2] - Q1[3]*Q2[3] - Q1[4]*Q2[4],
        Q1[1]*Q2[2] + Q1[2]*Q2[1] + Q1[3]*Q2[4] - Q1[4]*Q2[3],
        Q1[1]*Q2[3] - Q1[2]*Q2[4] + Q1[3]*Q2[1] + Q1[4]*Q2[2],
        Q1[1]*Q2[4] + Q1[2]*Q2[3] - Q1[3]*Q2[2] + Q1[4]*Q2[1])
end

Base.:/(Q1::Quaternion, Q2::Quaternion) = Q1 * inv(Q2)

Base.:*(Q::Quaternion, r::Real) = Quaternion(Q.coeffs * r)
Base.:*(r::Real, Q::Quaternion) = Quaternion(Q.coeffs * r)
Base.:/(Q::Quaternion, r::Real) = Quaternion(Q.coeffs / r)

Base.isreal(Q::Quaternion) = iszero(Q[2]) && iszero(Q[3]) && iszero(Q[4])
Base.isfinite(Q::Quaternion) = isfinite(Q[1]) && isfinite(Q[2]) && isfinite(Q[3]) && isfinite(Q[4])
Base.iszero(Q::Quaternion) = iszero(Q[1]) && iszero(Q[2]) && iszero(Q[3]) && iszero(Q[4])
Base.isone(Q::Quaternion) = isone(Q[1]) && iszero(Q[2]) && iszero(Q[3]) && iszero(Q[4])
Base.isnan(Q::Quaternion) = isnan(Q[1]) || isnan(Q[2]) || isnan(Q[3]) || isnan(Q[4])
Base.isinf(Q::Quaternion) = isinf(Q[1]) || isinf(Q[2]) || isinf(Q[3]) || isinf(Q[4])
Base.isinteger(Q::Quaternion) = isinteger(Q[1]) && isreal(Q)
