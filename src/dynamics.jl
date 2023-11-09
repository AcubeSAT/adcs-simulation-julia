qderiv(w, q) = 0.5 * (Quaternion(w...) * q)
wderiv(I, w, τ) = vec(I \ (τ - cross(w, (I * w))))

function wqderiv(I, w, τ, q)
    dw_dt = wderiv(I, w, τ)
    dq_dt = qderiv(w, q)
    return dw_dt, dq_dt
end

function rk4(I, w, τ, q, dt)
    k1_w, k1_q = wqderiv(I, w, τ, q)
    k2_w, k2_q = wqderiv(I, w .+ 0.5 .* dt .* k1_w, τ, q + 0.5 * dt * k1_q)
    k3_w, k3_q = wqderiv(I, w .+ 0.5 .* dt .* k2_w, τ, q + 0.5 * dt * k2_q)
    k4_w, k4_q = wqderiv(I, w .+ dt .* k3_w, τ, q + dt * k3_q)
    new_w = w .+ (dt / 6.0) .* (k1_w .+ 2 .* k2_w .+ 2 .* k3_w .+ k4_w)
    new_q = q + (dt / 6.0) * (k1_q + 2 * k2_q + 2 * k3_q + k4_q)
    new_q = normalize(new_q)
    return new_w, new_q
end

function rk4_filter(w, q, dt)
    k1_q = qderiv(w, q)
    k2_q = qderiv(w, q + 0.5 * dt * k1_q)
    k3_q = qderiv(w, q + 0.5 * dt * k2_q)
    k4_q = qderiv(w, q + dt * k3_q)
    new_q = normalize(q + (dt / 6.0) * (k1_q + 2 * k2_q + 2 * k3_q + k4_q))
    return new_q
end

function rk4_rw(J, w, τ, dt)
    k1_w = wderiv(J, w, τ)
    k2_w = wderiv(J, w + 0.5 * dt * k1_w, τ)
    k3_w = wderiv(J, w + 0.5 * dt * k2_w, τ)
    k4_w = wderiv(J, w + dt * k3_w, τ)
    return w + (dt / 6.0) * (k1_w + 2 * k2_w + 2 * k3_w + k4_w)
end

@concrete struct SimulationParams
    PD::PDController
    qtarget
    wtarget
    gr_model
    max_degree
    P
    dP
    msaturation
    sensors
    I
    dt
    m
end

@concrete struct SimulationContext
    state
    τw
    τsm
    τgravs
    τrmds
    magnetic_field_eci
    RW::ReactionWheel
end

# qtarget must be orbit2body otherwise I'll kick a hole in your fence
# TODO: what if saturation compensation is smaller than the cubesat w from control
function control_loop(Mode, SimParams::SimulationParams, SimContext::SimulationContext, PointingArgs::PointingArguments, r_ecef, epc, R_ecef_to_eci, mag_body, target_vectors, curindex)
    w, qeci2body = SimContext.state[curindex]
    qerr = emulate_estimation(SimParams.sensors, target_vectors, w)
    qestimated = qerr * mode_quaternion(Mode, PointingArgs)
    τ = calculate_torque(SimParams.PD, SimParams.qtarget, qestimated, w, SimParams.wtarget, qeci2body)
    τw, τsm, mtrue = decompose_torque(τ, mag_body, SimParams.msaturation)
    compensation = deadzone_compensation(SimContext.RW) + saturation_compensation(SimContext.RW)
    τw = clamp.(τw + compensation, -SimContext.RW.maxtorque, SimContext.RW.maxtorque)
    # rwfriction = stribeck(RW)
    @reset SimContext.RW.w = rk4_rw(SimContext.RW.J, SimContext.RW.w, -τw, SimParams.dt)
    τrmd = residual_dipole(SimParams.m, mag_body)
    G_ecef = gravity_gradient_tensor(SimParams.gr_model, r_ecef, epc, SimParams.max_degree, SimParams.P, SimParams.dP)
    R_ecef_to_body = to_rotation_matrix(qeci2body) * SMatrix{3,3}(R_ecef_to_eci)
    τgrav = gravity_torque(G_ecef, R_ecef_to_body, SimParams.I)
    wqeci2body = rk4(SimParams.I, w, τw + τsm + τrmd + τgrav, qeci2body, SimParams.dt)
    return wqeci2body, τw, τsm, τgrav, τrmd
end