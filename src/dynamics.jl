function derivatives(I, w, torque, q)
    # Calculate angular velocity derivatives in one line
    dw_dt = vec(I\(torque - w × (I * w)))

    # Calculate quaternion derivative
    dq_dt = 0.5 * (quat_mult(vcat(0,w),q))
    
    return dw_dt, dq_dt
end

# 4th Order Runge-Kutta (RK4) Method for both w and q
function rk4(I, w, torque, q, dt)
    k1_w, k1_q = derivatives(I, w, torque, q)
    k2_w, k2_q = derivatives(I, w .+ 0.5 .* dt .* k1_w, torque, q .+ 0.5 .* dt .* k1_q)
    k3_w, k3_q = derivatives(I, w .+ 0.5 .* dt .* k2_w, torque, q .+ 0.5 .* dt .* k2_q)
    k4_w, k4_q = derivatives(I, w .+ dt .* k3_w, torque, q .+ dt .* k3_q)

    new_w = w .+ (dt / 6.0) .* (k1_w .+ 2 .* k2_w .+ 2 .* k3_w .+ k4_w)
    new_q = q .+ (dt / 6.0) .* (k1_q .+ 2 .* k2_q .+ 2 .* k3_q .+ k4_q)
    
    # Normalize quaternion to ensure it stays a unit quaternion
    new_q = normalize(new_q)
    
    return new_w, new_q
end

function rk4_filter(I, w, torque, q, dt)
    _, k1_q = derivatives(I, w, torque, q)
    _, k2_q = derivatives(I, w, torque, q .+ 0.5 .* dt .* k1_q)
    _, k3_q = derivatives(I, w, torque, q .+ 0.5 .* dt .* k2_q)
    _, k4_q = derivatives(I, w, torque, q .+ dt .* k3_q)
    new_q = q .+ (dt / 6.0) .* (k1_q .+ 2 .* k2_q .+ 2 .* k3_q .+ k4_q)
    # Normalize quaternion to ensure it stays a unit quaternion
    new_q = normalize(new_q)
    return new_q
end