@concrete struct KalmanFilter
    transition_fun::Any
    transition_fun_jacobian::Any
    Q::Any
    measurement_fun::Any
    measurement_fun_jacobian::Any
    R::Any
    dt::Any
end

@concrete struct KFState
    q::Any
    bias::Any
end

function predict(state::KFState, P, KF::KalmanFilter, gyroscope_measurement)
    bias = state.bias
    F_k = KF.transition_fun_jacobian(gyroscope_measurement, bias)
    new_state = KF.transition_fun(state, gyroscope_measurement, KF.dt)
    Phi = exp(F_k * KF.dt)
    new_P = Phi * P * transpose(Phi) + KF.Q
    return (new_state, new_P)
end

function update(state::KFState,
    P,
    KF::KalmanFilter,
    groundtruth_measurements::Tuple,
    reference_vectors::Tuple)
    q = state.q
    mag_eci = reference_vectors[1]
    sun_eci = reference_vectors[2]
    H_k = KF.measurement_fun_jacobian(q, mag_eci, sun_eci)
    mag_body, sun_body = KF.measurement_fun(q, mag_eci, sun_eci)
    Kg = P * transpose(H_k) / (H_k * P * transpose(H_k) + KF.R)
    local_error_state = Kg * ([groundtruth_measurements[1]; groundtruth_measurements[2]] -
                         [mag_body; sun_body])
    local_error_quaternion = Quaternion([1; 0.5 * local_error_state[1:3]])
    new_q = normalize(q * local_error_quaternion)
    new_bias = state.bias + local_error_state[4:6]
    N_params = 6
    identity_matrix = 1.0 * I(N_params)
    new_P = (identity_matrix .- Kg * H_k) * P
    new_state = KFState(new_q, new_bias)
    return (new_state, new_P)
end

function transition_function(state::KFState, gyroscope_measurement, dt)
    q = state.q
    bias = state.bias
    w = gyroscope_measurement - bias
    q = rk4_filter(w, q, dt)
    return KFState(q, bias)
end

function transition_function_jacobian(gyroscope_measurement, bias)
    w = gyroscope_measurement - bias
    return vcat(hcat(-skew_symmetric(w), -1.0 * I(3)), zeros(3, 6))
end

function measurement_function(q, mag_eci, sun_eci)
    mag_body = rotvec(mag_eci, q)
    sun_body = rotvec(sun_eci, q)
    return (mag_body, sun_body)
end

function measurement_fun_jacobian(q, mag_eci, sun_eci)
    mag_body, sun_body = measurement_function(q, mag_eci, sun_eci)
    return [skew_symmetric(mag_body) zeros(3, 3); skew_symmetric(sun_body) zeros(3, 3)]
end
