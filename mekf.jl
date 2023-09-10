using LinearAlgebra
include("utils.jl")

mutable struct KalmanFilter
    # global_state::Vector{Float64}  # Global state vector  [1:4 quaternion 5:7 bias]
    # P::Matrix{Float64}  # State covariance matrix
    transition_fun::Function  # State transition model
    transition_fun_jacobian::Function
    Q::Matrix{Float64}  # Process noise covariance
    measurement_fun::Function  # Observation model
    measurement_fun_jacobian::Function
    R::Matrix{Float64}  # Observation noise covariance
    dt::Float64
end

function predict(state, P, kf::KalmanFilter, gyroscope_measurement)
    bias = state[5:7]
    q = state[1:4]
    F_k = kf.transition_fun_jacobian(gyroscope_measurement, bias)
    kf.transition_fun(q, gyroscope_measurement, bias, dt)
    new_state = kf.transition_fun(q, gyroscope_measurement, bias, kf.dt)
    # if new_state[1] < 0.0
    #     new_state[1:4] = -new_state[1:4]
    # end
    Phi = exp(F_k * kf.dt)
    new_P = Phi *  P * transpose(Phi) + kf.Q
    return (new_state, new_P)
end

function update(state, P, kf::KalmanFilter, groundtruth_measurements::Tuple, reference_vectors::Tuple)
    q = state[1:4]
    mag_eci = reference_vectors[1]
    sun_eci = reference_vectors[2]
    H_k = kf.measurement_fun_jacobian(state[1:4], mag_eci,sun_eci)
    mag_body, sun_body = kf.measurement_fun(q, mag_eci, sun_eci)
    Kg = P * transpose(H_k) / (H_k * P * transpose(H_k) + kf.R)
        
    local_error_state = Kg * ([groundtruth_measurements[1];groundtruth_measurements[2]] - [mag_body;sun_body])
    local_error_quaternion = [1;0.5 * local_error_state[1:3]]
    
    new_q = quat_mult(state[1:4],local_error_quaternion) 
    new_q = new_q/norm(new_q)
    new_bias = state[5:7] + local_error_state[4:6]            
    
    N_params = 6
    identity_matrix = 1.0 * I(N_params)
    
    new_P = (identity_matrix .- Kg * H_k) * P
    new_state = [new_q;new_bias]
    # if new_state[1] < 0.0
    #     new_state[1:4] = -new_state[1:4]
    # end
    return (new_state, new_P)
end

function transition_function(q, gyroscope_measurement, bias, dt)
    w = gyroscope_measurement - bias
    q = quat_mult(q,quaternion_exp([0;w]*dt))
    # torque = zeros(3)
    # dummy_I = I(3)
    # q = rk4_filter(I, w, torque, q, dt)
    next_state = vcat(q,bias)
end

function transition_function_jacobian(gyroscope_measurement, bias)
    w = gyroscope_measurement - bias
    F = vcat(hcat(-skew_symmetric(w),-1.0*I(3)), zeros(3,6))
end

function measurement_function(q, mag_eci, sun_eci)
    mag_body = rotate_vector_by_quaternion(mag_eci, q)
    sun_body = rotate_vector_by_quaternion(sun_eci, q)
    return (mag_body, sun_body)
end

function measurement_fun_jacobian(q, mag_eci, sun_eci)
    mag_body, sun_body = measurement_function(q, mag_eci, sun_eci)
    H = [skew_symmetric(mag_body) zeros(3,3);skew_symmetric(sun_body) zeros(3,3)]
end