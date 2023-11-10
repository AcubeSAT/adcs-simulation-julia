# Inspired by Michael Boyle and Base.Complex

# Signify standard arithmetic operations should be implemented on it
struct Quaternion{T<:Real} <: Number
    coeffs::SVector{4,T}
    Quaternion{T}(v::SVector{4,T}) where {T<:Real} = new{T}(v)
    # Always use SVector internally
    Quaternion{T}(v::A) where {T<:Real,A<:AbstractVector} = new{T}(SVector{4,T}(v))
end

# Construct without having to specify T
Quaternion(v::SVector{4,T}) where {T<:Real} = Quaternion{T}(v)
Quaternion(v::AbstractVector{T}) where {T<:Real} = Quaternion{T}(v)

# Explicitly state type, use SVector promote_rules for mixed types
Quaternion{T}(w, x, y, z) where {T<:Real} = Quaternion(SVector{4,T}(w, x, y, z))
Quaternion{T}(x, y, z) where {T<:Real} = Quaternion(SVector{4,T}(zero(x), x, y, z))
Quaternion{T}(w::Real) where {T<:Real} =
    Quaternion(SVector{4,T}(w, zero(w), zero(w), zero(w)))
# Rely on type inference, use SVector promote_rules for mixed types
function Quaternion(w, x, y, z)
    v = SVector{4}(w, x, y, z)
    return Quaternion{eltype(v)}(v)
end

function Quaternion(x, y, z)
    v = SVector{4}(zero(x), x, y, z)
    return Quaternion{eltype(v)}(v)
end

function Quaternion(w::Real)
    v = SVector{4}(w, zero(w), zero(w), zero(w))
    return Quaternion{eltype(v)}(v)
end

function Quaternion(axis::AbstractVector, angle::Real)
    h = angle / 2
    hsin = sin(h)
    q1 = cos(h)
    q2 = axis[1] * hsin
    q3 = axis[2] * hsin
    q4 = axis[3] * hsin
    return Quaternion(q1, q2, q3, q4)
end

# Type-preserving copy constructor
Quaternion{T}(Q::Quaternion{T}) where {T<:Real} = Quaternion(Q.coeffs...)
# Type conversion copy constructor
Quaternion{T}(Q::Quaternion{S}) where {T<:Real,S<:Real} =
    Quaternion(SVector{4,T}(Q.coeffs...))

const QuaternionF64 = Quaternion{Float64}
const QuaternionF32 = Quaternion{Float32}
const QuaternionF16 = Quaternion{Float16}

Base.zero(::Type{Quaternion{T}}) where {T<:Real} =
    Quaternion{T}(zero(T), zero(T), zero(T), zero(T))
Base.zero(Q::Quaternion{T}) where {T<:Real} = Base.zero(typeof(Q))

Base.one(::Type{Quaternion{T}}) where {T<:Real} =
    Quaternion{T}(one(T), zero(T), zero(T), zero(T))
Base.one(Q::Quaternion{T}) where {T<:Real} = Base.one(typeof(Q))

# Given a numeric type, return a Quaternion (not instance) specialized on T
(::Type{Quaternion})(::Type{T}) where {T<:Real} = Quaternion{T}
# Given a Quaternion specialized on T, return another Quaternion also specialized on T
(::Type{Quaternion})(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}

Base.eltype(::Type{Quaternion{T}}) where {T<:Real} = T
Base.promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T<:Real,S<:Real} =
    Quaternion{promote_type(T, S)}
Base.promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T<:Real,S<:Real} =
    Quaternion{promote_type(T, S)}

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
        Q1[1] * Q2[1] - Q1[2] * Q2[2] - Q1[3] * Q2[3] - Q1[4] * Q2[4],
        Q1[1] * Q2[2] + Q1[2] * Q2[1] + Q1[3] * Q2[4] - Q1[4] * Q2[3],
        Q1[1] * Q2[3] - Q1[2] * Q2[4] + Q1[3] * Q2[1] + Q1[4] * Q2[2],
        Q1[1] * Q2[4] + Q1[2] * Q2[3] - Q1[3] * Q2[2] + Q1[4] * Q2[1],
    )
end

Base.:/(Q1::Quaternion, Q2::Quaternion) = Q1 * inv(Q2)

Base.:*(Q::Quaternion, r::Real) = Quaternion(Q.coeffs * r)
Base.:*(r::Real, Q::Quaternion) = Quaternion(Q.coeffs * r)
Base.:/(Q::Quaternion, r::Real) = Quaternion(Q.coeffs / r)

Base.isreal(Q::Quaternion) = iszero(Q[2]) && iszero(Q[3]) && iszero(Q[4])
Base.isfinite(Q::Quaternion) =
    isfinite(Q[1]) && isfinite(Q[2]) && isfinite(Q[3]) && isfinite(Q[4])
Base.iszero(Q::Quaternion) = iszero(Q[1]) && iszero(Q[2]) && iszero(Q[3]) && iszero(Q[4])
Base.isone(Q::Quaternion) = isone(Q[1]) && iszero(Q[2]) && iszero(Q[3]) && iszero(Q[4])
Base.isnan(Q::Quaternion) = isnan(Q[1]) || isnan(Q[2]) || isnan(Q[3]) || isnan(Q[4])
Base.isinf(Q::Quaternion) = isinf(Q[1]) || isinf(Q[2]) || isinf(Q[3]) || isinf(Q[4])
Base.isinteger(Q::Quaternion) = isinteger(Q[1]) && isreal(Q)

function Random.rand(::Type{Quaternion})
    elements = normalize(rand(4))
    return Quaternion(elements)
end

function rotvec(v::A, Q::Quaternion) where {A<:AbstractVector}
    Q = normalize(Q)
    return vec(Q * Quaternion(v[1], v[2], v[3]) * conj(Q))
end

# From https://github.com/moble/Quaternionic.jl
dominant_eigenvector(M::Symmetric{T,SMatrix{4,4,T,16}}) where {T} = eigen(M).vectors[:, 4]
dominant_eigenvector(M::Symmetric{Float64,SMatrix{4,4,Float64,16}}) =
    eigen(M, 4:4).vectors[:, 1]
dominant_eigenvector(M::Symmetric{Float32,SMatrix{4,4,Float32,16}}) =
    eigen(M, 4:4).vectors[:, 1]

"""
    from_rotation_matrix(A)

Convert 3x3 rotation matrix to quaternion.

Assuming the 3x3 matrix `A` rotates a vector `v` according to

    v' = A * v,

we can also express this rotation in terms of a quaternion `R` such that

    v' = R * v * R⁻¹.

This function returns that quaternion, using Bar-Itzhack's algorithm (version 3) to allow
for non-orthogonal matrices.  [J. Guidance, Vol. 23, No. 6, p.
1085](http://dx.doi.org/10.2514/2.4654)

!!! note
    If you want to use this function for matrices with elements of types other than
    `Float64` or `Float32`, you will need to (install and) import `GenericLinearAlgebra`
    first.  The reason is that this function computes the eigen-decomposition of `A`, which
    is only available for more generic float types via that package.  Note that you will
    want at least version 0.3.11 of `GenericLinearAlgebra` because previous versions had a
    bug.

"""
function from_rotation_matrix(A::AbstractMatrix)
    @assert size(A) == (3, 3)
    @inbounds begin
        # Compute 3K₃ according to Eq. (2) of Bar-Itzhack.  We will just be looking for the
        # eigenvector with the largest eigenvalue, so scaling by a strictly positive number
        # (3, in this case) won't change that.
        K = Symmetric(
            @SMatrix[
                A[1, 1]-A[2, 2]-A[3, 3] A[2, 1]+A[1, 2] A[3, 1]+A[1, 3] A[2, 3]-A[3, 2]
                0 A[2, 2]-A[1, 1]-A[3, 3] A[3, 2]+A[2, 3] A[3, 1]-A[1, 3]
                0 0 A[3, 3]-A[1, 1]-A[2, 2] A[1, 2]-A[2, 1]
                0 0 0 A[1, 1]+A[2, 2]+A[3, 3]
            ]
        )

        # Compute the *dominant* eigenvector (the one with the largest eigenvalue)
        de = dominant_eigenvector(K)

        # Convert it into a quaternion
        R = Quaternion(de[4], -de[1], -de[2], -de[3])
    end
    R
end

"""
    to_rotation_matrix(q)

Convert quaternion to 3x3 rotation matrix.

Assuming the quaternion `R` rotates a vector `v` according to

    v' = R * v * R⁻¹,

we can also express this rotation in terms of a 3x3 matrix `ℛ` such that

    v' = ℛ * v.

This function returns that matrix.

"""
function to_rotation_matrix(q::Quaternion)
    n = inv(abs2(q))
    @SMatrix [
        1-2*(q[3]^2+q[4]^2)*n 2*(q[2]*q[3]-q[4]*q[1])*n 2*(q[2]*q[4]+q[3]*q[1])*n
        2*(q[2]*q[3]+q[4]*q[1])*n 1-2*(q[2]^2+q[4]^2)*n 2*(q[3]*q[4]-q[2]*q[1])*n
        2*(q[2]*q[4]-q[3]*q[1])*n 2*(q[3]*q[4]+q[2]*q[1])*n 1-2*(q[2]^2+q[3]^2)*n
    ]
end

function to_euler_angles(q::Quaternion)
    a0 = 2acos(√((q[1]^2 + q[4]^2) / abs2(q)))
    a1 = atan(q[4], q[1])
    a2 = atan(-q[2], q[3])
    @SVector [a1 + a2, a0, a1 - a2]
end

function euler_to_quaternion(roll::Real, pitch::Real, yaw::Real)
    half_roll = roll / 2
    half_pitch = pitch / 2
    half_yaw = yaw / 2

    cr, sr = cos(half_roll), sin(half_roll)
    cp, sp = cos(half_pitch), sin(half_pitch)
    cy, sy = cos(half_yaw), sin(half_yaw)

    w = cr * cp * cy + sr * sp * sy
    x = sr * cp * cy - cr * sp * sy
    y = cr * sp * cy + sr * cp * sy
    z = cr * cp * sy - sr * sp * cy

    Quaternion(w, x, y, z)
end    

function Base.show(io::IO, Q::Quaternion{T}) where {T}
    print(io, "{")
    show(io, T)
    print(io, "}")
    print(io, " q1: ")
    show(io, Q[1])
    print(io, ", q2: ")
    show(io, Q[2])
    print(io, ", q3: ")
    show(io, Q[3])
    print(io, ", q4: ")
    show(io, Q[4])
end
