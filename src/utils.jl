function quat_mult(q, p)
    return [
        q[1] * p[1] - q[2] * p[2] - q[3] * p[3] - q[4] * p[4],
        q[1] * p[2] + q[2] * p[1] + q[3] * p[4] - q[4] * p[3],
        q[1] * p[3] - q[2] * p[4] + q[3] * p[1] + q[4] * p[2],
        q[1] * p[4] + q[2] * p[3] - q[3] * p[2] + q[4] * p[1],
    ]
end

function quaternion_conjugate(q::Vector{Float64})
    w, x, y, z = q
    return [w, -x, -y, -z]
end

function rotate_vector_by_quaternion(v::AbstractVector, q::Vector{Float64})
    v_quaternion = [0.0; v]

    v_rotated_quaternion = quat_mult(quaternion_conjugate(q), quat_mult(v_quaternion, q))

    v_rotated = v_rotated_quaternion[2:4]

    return v_rotated
end

function skew_symmetric(v::AbstractArray)
    x, y, z = v
    S = [0.0 -z y;
        z 0.0 -x;
        -y x 0.0]
    return S
end

function package_weights(x)
    Q = diagm(exp.([x[1], x[1], x[1], x[2], x[2], x[2]]))
    R = diagm(exp.([x[3], x[3], x[3], x[4], x[4], x[4]]))
    return (Q, R)
end
