function calculate_orbit(JD, n_orbits, dt)
    epc0 = Epoch(jd_to_caldate(JD)...)
    oe0 = [R_EARTH + 500e3, 0.01, 75.0, 45.0, 30.0, 0.0]
    eci0 = sOSCtoCART(oe0, use_degrees = true)
    T = orbit_period(oe0[1])
    epcf = epc0 + n_orbits * T
    orb = EarthInertialState(epc0, eci0, dt = dt,
        mass = 100.0, n_grav = 20, m_grav = 20,
        drag = true, srp = true,
        moon = true, sun = true,
        relativity = false)

    t, epc, eci = sim!(orb, epcf)
    return t, epc, eci
end

function run_groundtruth_simulation(params)
    torque = zeros(3, 1)
    w = 0.0035 * ones(3)
    q = Quaternion([1.0, 0, 0, 0])
    JD = 2459921.0
    n_orbits = 2
    _, epc, eci = calculate_orbit(JD, n_orbits)
    r_eci = eci[1:3, :]
    r_eci = [r_eci[:, i] for i in 1:size(r_eci, 2)]
    rotation_eci2ecef = rECItoECEF.(epc)
    r_ecef = [rotation_eci2ecef[i] * r_eci[i] for i in 1:size(rotation_eci2ecef, 1)]
    rotation_ecef2eci = rECEFtoECI.(epc)
    mag_ecef = geomagnetic_dipole_field.(r_ecef)
    mag_ecef = [mag_ecef[i] ./ norm(mag_ecef[i]) for i in 1:size(mag_ecef, 1)]
    mag_eci = [rotation_ecef2eci[i] * mag_ecef[i] for i in 1:size(rotation_ecef2eci, 1)]
    sun_eci = sun_position.(epc)
    sun_eci = [sun_eci[i] ./ norm(sun_eci[i]) for i in 1:size(sun_eci, 1)]
    w_history = Array{Float64}(undef, 3, length(epc))
    q_history = Vector{Quaternion}(undef, length(epc))
    bias_history = Array{Float64}(undef, 3, length(epc))
    for i in 1:length(epc)
        w, q = rk4(params.inertia_matrix, w, torque, q, params.dt)
        w_history[:, i] = w
        q_history[i] = q
    end
    bias = 0.01 * ones(3)
    mag_noisy_history = Array{Float64}(undef, 3, length(epc))
    sun_noisy_history = Array{Float64}(undef, 3, length(epc))
    gyro_noisy_history = Array{Float64}(undef, 3, length(epc))

    for i in 1:length(epc)
        mag_noisy_history[:, i], sun_noisy_history[:, i], gyro_noisy_history[:, i], bias = get_noisy_measurements(q_history[i],
            w_history[:, i],
            bias,
            mag_eci[i],
            sun_eci[i],
            params)

        bias_history[:, i] = bias
    end
    groundtruth_state_history = (q_history,
        w_history,
        bias_history,
        mag_eci,
        sun_eci,
        mag_noisy_history,
        sun_noisy_history,
        gyro_noisy_history)

    return groundtruth_state_history
end

function run_filter_simulation(tunable_params,
    params,
    mag_noisy,
    sun_noisy,
    mag_eci,
    sun_eci,
    gyroscope_measurement)
    kf = KalmanFilter(transition_function,
        transition_function_jacobian,
        tunable_params[1],
        measurement_function,
        measurement_fun_jacobian,
        tunable_params[2],
        params.dt)

    state = KFState(Quaternion(1.0, 0.0, 0.0, 0.0), zeros(3))
    P = 1.0 * Matrix{Float64}(I, 6, 6)
    N = size(mag_noisy, 2)
    state_estimation_array = Vector{KFState}(undef, N)
    for i in 1:size(mag_noisy)[2]
        state, P = update(state,
            P,
            kf,
            (mag_noisy[:, i], sun_noisy[:, i]),
            (mag_eci[i], sun_eci[i]))

        state, P = predict(state, P, kf, gyroscope_measurement[:, i])
        state_estimation_array[i] = state
    end
    return state_estimation_array
end

# FOV should be in radians, half of the sensor FOV
function in_fov(tpos, fov)
    LinearAlgebra.normalize!(tpos)
    return acos(dot(tpos, [0 0 1])) <= fov
end

function in_fov(tpos, vfov, hfov)
    LinearAlgebra.normalize!(tpos)
    θv = abs(acos(tpos[3]))
    θh = abs(atan(tpos[2], tpos[1]))
    return θv <= vfov && θh <= hfov
end

function available(NS::NadirSensor, target_vector, w)
    return w <= NS.maximum_rate && in_fov(target_vector, NS.vfov, NS.hfov)
end

function available(ST::StarTracker, target_vector, w)
    return w <= ST.maximum_rate && !in_fov(target_vector, ST.fov)
end

function available(SN::SunSensor, target_vector, w)
    return w <= SN.maximum_rate && in_fov(target_vector, SN.fov)
end

available(::AbstractSensor) = error("available is not defined in the abstract type")
function create_err_q(S::AbstractSensor)
    δθx = randn() * S.σ_cross_sight
    δθy = randn() * S.σ_cross_sight
    δθz = randn() * S.σ_roll
    return LinearAlgebra.normalize(Quaternion(1.0, δθx / 2, δθy / 2, δθz / 2))
end

function estimateq(sensors, target_vectors, w)
    err_qs = Quaternion[]
    abs_weights = typeof(sensors[1].abs_weight)[]
    for (sensor, tvec) in zip(sensors, target_vectors)
        if available(sensor, tvec, w)
            push!(err_qs, create_err_q(sensor))
            push!(abs_weights, sensor.abs_weight)
        end
    end
    isempty(err_qs) && return nothing
    rel_weights = abs_weights ./ sum(abs_weights)
    return normalize(sum(err_qs .* rel_weights))
end
