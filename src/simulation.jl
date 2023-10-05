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
    _, epc, eci = calculate_orbit(JD, n_orbits, 0.001)
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

#TODO rw_w
function rotational_dynamics(qeci2body, w, pointing_mode, PD, t, epc::Vector{Epoch}, r_eci, r_ecef, sun_eci, mag_eci, Qeci2orbit, R_ecef_to_eci, dt, qtarget, model, max_degree, P, dP)
    rw_w = 94.247779 * ones(3)
    sensors = (NadirSensor(), StarTracker(), SunSensor())
    RW = ReactionWheel(J=I(3), w=rw_w, saturationα=1, deadzoneα=1, maxtorque=0.001)
    iters = length(t)
    state_history = Vector{Tuple{typeof(w), typeof(qeci2body)}}(undef, iters)
    τw = Vector{Vector{Float64}}(undef, iters)
    τsm = Vector{Vector{Float64}}(undef, iters)
    res = Vector{Bool}(undef, iters)
    τgravs = Vector{Vector{Float64}}(undef, iters)
    τrmds = Vector{Vector{Float64}}(undef, iters)

    t = epoch_to_datetime(epc)
    for i in 1:iters
        nadir_body = Vector(-rotvec(normalize(r_eci[i]), qeci2body))
        sun_body = Vector(rotvec(sun_eci[i], qeci2body))
        mag_body = Vector(rotvec(mag_eci[i], qeci2body))
        qeci2orbit = Qeci2orbit[i]
        PA = PointingArguments(sun_body, nadir_body, qeci2body, qeci2orbit)
        wqeci2body, rτw, rτsm, rres, RW, τgrav, τrmd = control_loop(pointing_mode, PA, PD, qeci2body, qeci2orbit, qtarget, zeros(3), mag_body, 0.66, sensors, (nadir_body, sun_body, sun_body), w, RW, diagm([0.167,0.142,0.142]), model, r_ecef[i], t[i], max_degree, P, dP, R_ecef_to_eci[i], dt)
        w, qeci2body = wqeci2body
        state_history[i] = (w, qeci2body)
        τw[i] = rτw
        τsm[i] = rτsm
        res[i] = rres
        τgravs[i] = τgrav
        τrmds[i] = τrmd
    end
    return state_history, τw, τsm, τgravs, τrmds
end

function generate_orbit_data(jd, norbits, dt)
    t, epc, eci = calculate_orbit(jd, norbits, dt)
    r_eci = eci[1:3, :]
    r_eci = [r_eci[:, i] for i in 1:size(r_eci, 2)]
    v_eci = eci[4:6, :]
    v_eci = [v_eci[:, i] for i in 1:size(v_eci, 2)]
    rotation_eci2ecef = rECItoECEF.(epc)
    r_ecef = [rotation_eci2ecef[i] * r_eci[i] for i in 1:size(rotation_eci2ecef, 1)]
    rotation_ecef2eci = rECEFtoECI.(epc)
    mag_ecef = 1e-9 * geomagnetic_dipole_field.(r_ecef)
    mag_ecef = [mag_ecef[i] for i in 1:size(mag_ecef, 1)]
    mag_eci = [rotation_ecef2eci[i] * mag_ecef[i] for i in 1:size(rotation_ecef2eci, 1)]
    sun_eci = sun_position.(epc)
    sun_eci = [normalize(sun_eci[i]) for i in 1:size(sun_eci, 1)]
    T = eci2orbit.(r_eci, v_eci)
    qeci2orbit = from_rotation_matrix.(T)
    return t, epc, r_eci, r_ecef, sun_eci, mag_eci, qeci2orbit, rotation_ecef2eci
end

function eci2orbit(r_eci, v_eci)
    z = -normalize(r_eci)
    y = normalize(cross(r_eci, v_eci))
    x = cross(y, z)
    T = hcat(x, y, z)
    return T'
end
