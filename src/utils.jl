function skew_symmetric(v::AbstractArray)
    x, y, z = v
    S = [0.0 -z y
        z 0.0 -x
        -y x 0.0]
    return S
end

function package_weights(x)
    Q = diagm(exp.([x[1], x[1], x[1], x[2], x[2], x[2]]))
    R = diagm(exp.([x[3], x[3], x[3], x[4], x[4], x[4]]))
    return (Q, R)
end

_struct_to_vec(S) = [getfield(S, field) for field in fieldnames(typeof(S))]
