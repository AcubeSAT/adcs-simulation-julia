using Statistics

mse(ŷ, y; agg=mean) = agg(abs2.(ŷ .- y))
