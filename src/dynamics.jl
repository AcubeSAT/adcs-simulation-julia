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

# qtarget must be orbit2body otherwise I'll kick a hole in your fence
# TODO: what if saturation compensation is smaller than the cubesat w from control
function control_loop(sun_tracking, PD, qeci2body, qeci2orbit, qtarget, wtarget, b, msaturation, sensors, target_vectors, w, RW::ReactionWheel, I, model, r_ecef, epc, max_degree, P, dP, R_ecef_to_eci, dt)
    res, qerr = emulate_estimation(sensors, target_vectors, w)
    qorbit2body = qeci2body * conj(qeci2orbit)
    qbody2sun = align_frame_with_vector(target_vectors[2], target_vectors[1], [0,0,-1],[0,1,0])
    if sun_tracking
        qestimated = qerr * qbody2sun
    else
        qestimated = qerr * qorbit2body
    end
    τ = calculate_torque(PD, qtarget, qestimated, w, wtarget, qeci2body)
    τw, τsm, mtrue = decompose_torque(τ, b, msaturation)
    compensation = deadzone_compensation(RW) + saturation_compensation(RW)
    τw = clamp.(τw + compensation, -RW.maxtorque, RW.maxtorque)
    # rwfriction = stribeck(RW)
    @reset RW.w = rk4_rw(RW.J, RW.w, -τw, dt)
    τrmd = residual_dipole([-0.1235, 0.2469, -0.2273], b)
    G_ecef = gravity_gradient_tensor(model, r_ecef, epc, max_degree, P, dP)
    R_ecef_to_body = to_rotation_matrix(qeci2body) * R_ecef_to_eci
    τgravity = gravity_torque(G_ecef, R_ecef_to_body, I)
    return rk4(I, w, τw + τsm + τrmd + τgravity, qeci2body, dt), τw, τsm, res, RW, τgravity, τrmd # TODO: update cubesat state with τw + τsm + disturbances
end
