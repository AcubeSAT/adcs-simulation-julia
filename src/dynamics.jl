qderiv(w, q) = 0.5 * (Quaternion([w; 0.0]) * q)

function wqderiv(I, w, torque, q)
    dw_dt = vec(I \ (torque - w × (I * w)))
    dq_dt = qderiv(w, q)
    return dw_dt, dq_dt
end

function rk4(I, w, torque, q, dt)
    k1_w, k1_q = wqderiv(I, w, torque, q)
    k2_w, k2_q = wqderiv(I, w .+ 0.5 .* dt .* k1_w, torque, q + 0.5 * dt * k1_q)
    k3_w, k3_q = wqderiv(I, w .+ 0.5 .* dt .* k2_w, torque, q + 0.5 * dt * k2_q)
    k4_w, k4_q = wqderiv(I, w .+ dt .* k3_w, torque, q + dt * k3_q)
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
