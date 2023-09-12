@concrete struct Quaternion
    q1
    q2
    q3
    q4
end

Quaternion(xs) = Quaternion(xs...)

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

function Base.:*(Q::Quaternion, r::Real)
    return Quaternion(Q.q1 * r, Q.q2 * r, Q.q3 * r, Q.q4 * r)
end

function Base.:*(r::Real, Q::Quaternion)
    return Q * r
end

function Base.:+(Q1::Quaternion, Q2::Quaternion)
    return Quaternion(Q1.q1 + Q2.q1, Q1.q2 + Q2.q2, Q1.q3 + Q2.q3, Q1.q4 + Q2.q4)
end

function Base.conj(Q::Quaternion)
    return Quaternion(Q.q1, -Q.q2, -Q.q3, -Q.q4)
end

function LinearAlgebra.norm(Q::Quaternion)
    return sqrt(Q.q1^2 + Q.q2^2 + Q.q3^2 + Q.q4^2)
end

function LinearAlgebra.normalize(Q::Quaternion)
    qnorm = norm(Q)
    return Quaternion(Q.q1 / qnorm, Q.q2 / qnorm, Q.q3 / qnorm, Q.q4 / qnorm)
end

scalar(Q::Quaternion) = Q.q1
vector(Q::Quaternion) = [Q.q2, Q.q3, Q.q4]

function rotvec(v::Vector{Float64}, Q::Quaternion)
    Q = normalize(Q)
    return vector(conj(Q) * Quaternion([0; v]) * Q)
end
