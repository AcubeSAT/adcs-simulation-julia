mse(ŷ, y; agg=mean) = agg(abs2.(ŷ .- y))

function qloss(q̂, q)
    relq = q * conj(q̂)
    return norm(vec(relq))
end