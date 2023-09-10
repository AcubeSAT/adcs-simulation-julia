using LinearAlgebra
include("utils.jl")

function get_noisy_measurements(q, gt_angular_velocity, bias, mag_ref, sun_ref, params)
    old_bias = bias
    bias = bias .+ params.sigma_u .* sqrt(params.dt) .* randn(3)
    total_noise =  0.5 .* (old_bias .+ bias) .+ sqrt((params.sigma_v^2 ./ params.dt .+ params.sigma_u^2 .* params.dt ./ 12)) .* randn(3)
    gyroscope_measurement = gt_angular_velocity + total_noise
    mag_noisy = rotate_vector_by_quaternion(mag_ref, q) + params.mag_noise * randn(3)
    sun_noisy = rotate_vector_by_quaternion(sun_ref, q) + params.sun_noise * randn(3)
    return (mag_noisy, sun_noisy, gyroscope_measurement, bias)
end