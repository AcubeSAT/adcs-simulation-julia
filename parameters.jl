using LinearAlgebra

mutable struct parameter_struct
    inertia_matrix::Matrix{Float64}
    sigma_u::Float64
    sigma_v::Float64
    mag_noise::Float64
    sun_noise::Float64
    dt::Float64
end