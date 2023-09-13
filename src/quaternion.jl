struct Quaternion{T}
    q1::T
    q2::T
    q3::T
    q4::T
end

function Base.show(io::IO, Q::Quaternion)
    print(io, "q1: ")
    show(io, Q.q1)
    print(io, ", q2: ")
    show(io, Q.q2)
    print(io, ", q3: ")
    show(io, Q.q3)
    print(io, ", q4: ")
    show(io, Q.q4)
end

Quaternion(xs::Vector) = Quaternion(xs...)
function Base.convert(::Type{Quaternion{T}}, x::T) where {T}
    return Quaternion(x, x, x, x)
end

function Base.convert(::Type{T}, Q::Quaternion) where {T}
    return Quaternion{T}(Q.q1, Q.q2, Q.q3, Q.q4)
end

Quaternion(x::T) where {T} = convert(Quaternion{T}, x)
Base.convert(::Type{Quaternion{T}}, x::Quaternion{T}) where {T} = x
Base.eltype(::Type{Quaternion{T}}) where {T} = T
Base.zero(Q::Quaternion{T}) where {T} = Quaternion(zero(eltype(Q)))
Base.length(::Quaternion) = 4
function Base.getindex(Q::Quaternion, i::Int)
    1 <= i <= 4 || throw(BoundsError(Q, i))
    i == 1 && return Q.q1
    i == 2 && return Q.q2
    i == 3 && return Q.q3
    i == 4 && return Q.q4
end

function Base.:*(Q1::Quaternion, Q2::Quaternion)
    return Quaternion(
        Q1.q1 * Q2.q1 - Q1.q2 * Q2.q2 - Q1.q3 * Q2.q3 - Q1.q4 * Q2.q4,
        Q1.q1 * Q2.q2 + Q1.q2 * Q2.q1 + Q1.q3 * Q2.q4 - Q1.q4 * Q2.q3,
        Q1.q1 * Q2.q3 - Q1.q2 * Q2.q4 + Q1.q3 * Q2.q1 + Q1.q4 * Q2.q2,
        Q1.q1 * Q2.q4 + Q1.q2 * Q2.q3 - Q1.q3 * Q2.q2 + Q1.q4 * Q2.q1
    )
end

Base.:*(Q::Quaternion, n::Number) = Quaternion(Q.q1 * n, Q.q2 * n, Q.q3 * n, Q.q4 * n)
Base.:*(n::Number, Q::Quaternion) = Q * n
Base.:/(Q::Quaternion, n::Number) = Quaternion(Q.q1 / n, Q.q2 / n, Q.q3 / n, Q.q4 / n)
function Base.:+(Q1::Quaternion, Q2::Quaternion)
    return Quaternion(Q1.q1 + Q2.q1, Q1.q2 + Q2.q2, Q1.q3 + Q2.q3, Q1.q4 + Q2.q4)
end

Base.:-(Q::Quaternion) = Quaternion(-Q.q1, -Q.q2, -Q.q3, -Q.q4)
function Base.:-(Q1::Quaternion, Q2::Quaternion)
    return Quaternion(Q1.q1 - Q2.q1, Q1.q2 - Q2.q2, Q1.q3 - Q2.q3, Q1.q4 - Q2.q4)
end

function Base.:(==)(Q1::Quaternion, Q2::Quaternion)
    return Q1.q1 == Q2.q1 && Q1.q2 == Q2.q2 && Q1.q3 == Q2.q3 && Q1.q4 == Q2.q4
end

Base.conj(Q::Quaternion) = Quaternion(Q.q1, -Q.q2, -Q.q3, -Q.q4)
LinearAlgebra.norm(Q::Quaternion) = sqrt(Q.q1^2 + Q.q2^2 + Q.q3^2 + Q.q4^2)
function LinearAlgebra.normalize(Q::Quaternion)
    qnorm = norm(Q)
    return Quaternion(Q.q1 / qnorm, Q.q2 / qnorm, Q.q3 / qnorm, Q.q4 / qnorm)
end

Base.:/(Q1::Quaternion, Q2::Quaternion) = Q1 * conj(Q2) / norm(Q2)^2
Base.inv(Q::Quaternion) = conj(Q) / norm(Q)^2
scalar(Q::Quaternion) = Q.q1
vector(Q::Quaternion) = [Q.q2, Q.q3, Q.q4]
function rotvec(v::Vector, Q::Quaternion)
    Q = normalize(Q)
    return vector(conj(Q) * Quaternion([0; v]) * Q)
end

function quat_to_euler_deg(Q::Quaternion)
    roll = atan(2 * (Q.q1 * Q.q2 + Q.q3 * Q.q4), 1 - 2 * (Q.q2^2 + Q.q3^2))
    pitch = asin(2 * (Q.q1 * Q.q3 - Q.q4 * Q.q2))
    yaw = atan(2 * (Q.q1 * Q.q4 + Q.q2 * Q.q3), 1 - 2 * (Q.q3^2 + Q.q4^2))
    return [(roll * 180 / π), (pitch * 180 / π), (yaw * 180 / π)]
end