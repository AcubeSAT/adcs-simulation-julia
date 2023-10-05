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

function split_into_parts(vec::Vector, n::Int)
    len = length(vec)
    size_per_part = div(len, n)

    parts = Vector{typeof(vec)}(undef, n)
    start_idx = 1

    for i in 1:n-1
        parts[i] = vec[start_idx:start_idx+size_per_part-1]
        start_idx += size_per_part
    end
    parts[n] = vec[start_idx:end] # last part takes the remainder

    return parts
end

function subvector(t::Tuple, start_idx, end_idx=:end)
    end_idx = end_idx == :end ? length(first(t)) : end_idx
    return Tuple(v[start_idx:end_idx] for v in t)
end
