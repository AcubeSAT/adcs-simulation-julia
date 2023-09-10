include("dynamics.jl")
include("utils.jl")
include("parameters.jl")
include("mekf.jl")
include("noise_models.jl")

function calculate_orbit(JD, n_orbits)
    # Declare simulation initial Epoch
    epc0 = Epoch(jd_to_caldate(JD)...)
    # Declare initial state in terms of osculating orbital elements
    oe0 = [R_EARTH + 500e3, 0.01, 75.0, 45.0, 30.0, 0.0]

    # Convert osculating elements to Cartesian state
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    # Set the propagation end time to one orbit period after the start
    T = orbit_period(oe0[1])
    epcf = epc0 + n_orbits * T

    # Initialize State Vector
    orb = EarthInertialState(epc0, eci0, dt=1.0,
        mass=100.0, n_grav=20, m_grav=20,
        drag=true, srp=true,
        moon=true, sun=true,
        relativity=false
    )

    # Propagate the orbit
    t, epc, eci = sim!(orb, epcf)
    return t, epc, eci
end


function run_groundtruth_simulation(params)
    torque = zeros(3, 1)
    w = 0.0035 * ones(3)
    q = [1.0, 0, 0, 0]
    JD = 2459921.0
    n_orbits = 2
    t, epc, eci = calculate_orbit(JD, n_orbits)

    r_eci = eci[1:3, :]
    r_eci = [r_eci[:, i] for i in 1:size(r_eci, 2)]
    # v_eci = eci[4:6,:]

    rotation_eci2ecef = rECItoECEF.(epc)
    r_ecef = [rotation_eci2ecef[i] * r_eci[i] for i in 1:size(rotation_eci2ecef, 1)]
    rotation_ecef2eci = rECEFtoECI.(epc)

    mag_ecef = geomagnetic_dipole_field.(r_ecef)
    mag_ecef = [mag_ecef[i] ./ norm(mag_ecef[i]) for i in 1:size(mag_ecef, 1)]
    mag_eci = [rotation_ecef2eci[i] * mag_ecef[i] for i in 1:size(rotation_ecef2eci, 1)]

    sun_eci = sun_position.(epc)
    sun_eci = [sun_eci[i] ./ norm(sun_eci[i]) for i in 1:size(sun_eci, 1)]

    w_history = Array{Float64}(undef, 3, length(epc))
    q_history = Array{Float64}(undef, 4, length(epc))
    bias_history = Array{Float64}(undef, 3, length(epc))

    # Dynamics

    for i in 1:length(epc)
        w, q = rk4(params.inertia_matrix, w, torque, q, params.dt)
        # dw_dt, dq_dt =  derivatives(I, w, torque, q)
        # q = q + dq_dt * params.dt
        # q = q/norm(q)
        # w = w + dw_dt * params.dt
        w_history[:, i] = w
        q_history[:, i] = q
    end

    bias = 0.001 * ones(3)
    mag_noisy_history = Array{Float64}(undef, 3, length(epc))
    sun_noisy_history = Array{Float64}(undef, 3, length(epc))
    gyro_noisy_history = Array{Float64}(undef, 3, length(epc))

    for i in 1:length(epc)
        mag_noisy_history[:, i], sun_noisy_history[:, i], gyro_noisy_history[:, i], bias = get_noisy_measurements(q_history[:, i], w_history[:, i], bias, mag_eci[i], sun_eci[i], params)
        bias_history[:, i] = bias
    end
    groundtruth_state_history = (q_history, w_history, bias_history, mag_eci, sun_eci, mag_noisy_history, sun_noisy_history, gyro_noisy_history)

    return groundtruth_state_history
end

function run_filter_simulation(tunable_params, params, mag_noisy, sun_noisy, mag_eci, sun_eci, gyroscope_measurement)

    kf = KalmanFilter(
        transition_function,             # F (state transition model)
        transition_function_jacobian,
        tunable_params[1],        # Q (process noise covariance)
        measurement_function, # H (observation model)
        measurement_fun_jacobian,
        tunable_params[2],        # R (observation noise covariance)
        params.dt             # I (identity matrix)
    )
    state = [1.0; 0.0; 0.0; 0; 0; 0; 0]
    P = 1.0 * Matrix{Float64}(I, 6, 6)

    N = size(mag_noisy, 2)
    state_estimation_array = Matrix{Float64}(undef, 7, N) # pre-allocate

    # Main loop
    for i in 1:size(mag_noisy_history)[2]
        state, P = update(state, P, kf, (mag_noisy[:, i], sun_noisy[:, i]), (mag_eci[i], sun_eci[i]))
        state, P = predict(state, P, kf, gyroscope_measurement[:, i])
        state_estimation_array[:, i] = state
    end
    return state_estimation_array
end

