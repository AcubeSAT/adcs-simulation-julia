qderiv(w, q) = 0.5 * (Quaternion([w; 0.0]) * q)
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

# TODO: what if saturation compensation is smaller than the cubesat w from control
# function control_loop(sensors, target_vectors, w, RW::ReactionWheel, I, dt = 0.001)
#     q = estimateq(sensors, target_vectors, w)
#     τ = calculate_torque(PDController(0.1, 0.001), qtarget, q, w, wtarget)
#     τw, τsm, mtrue = decompose_torque(τ, b, msaturation)
#     compensation = deadzone_compensation(RW) + saturation_compensation(RW)
#     τw = clamp(τw + compensation, -RW.maxtorque, RW.maxtorque)
#     rwfriction = stribeck(RW)
#     @reset RW.w = rk4_rw(RW.J, RW.w, -τw + rwfriction, dt)
#     return rk4(I, w, τw - rwfriction, q, dt) # TODO: update cubesat state with τw + τsm + disturbances
# end
