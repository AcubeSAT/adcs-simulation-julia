function plot_histories(groundtruth_state_history, state_estimation_history)
    p = plot(layout = (7, 1), size = (600, 1200))

    for i in 1:7
        plot!(p[i], groundtruth_state_history[i, :], label = "Ground Truth", xlabel = "Timestep", ylabel = "Value", linewidth = 2)
        plot!(p[i], state_estimation_history[i, :], label = "Estimated", xlabel = "Timestep", ylabel = "Value", linewidth = 2)
    end

    display(current())
end

function quat_to_euler_deg(q::Array{Float64, 1})
    qw, qx, qy, qz = q

    roll = atan(2 * (qw * qx + qy * qz), 1 - 2 * (qx^2 + qy^2))
    pitch = asin(2 * (qw * qy - qz * qx))
    yaw = atan(2 * (qw * qz + qx * qy), 1 - 2 * (qy^2 + qz^2))
    
    return [(roll * 180 / π), (pitch * 180 / π), (yaw * 180 / π)]
end

function plot_difference(gt_target::Array{Float64, 2}, state_estimation_history::Array{Float64, 2})
    if size(gt_target) != size(state_estimation_history)
        println("Error: Dimensions of gt_target and state_estimation_history must match.")
        return
    end
    
    nrows, ncols = size(gt_target)
    
    difference_array = zeros(Float64, nrows, ncols)
    
    for t in 1:ncols
        q_gt = gt_target[1:4, t]
        q_est = state_estimation_history[1:4, t]
        
        q_gt ./= norm(q_gt)
        q_est ./= norm(q_est)
        
        q_rel = quat_mult(q_est, quaternion_conjugate(q_gt))
        
        euler_angles_deg = quat_to_euler_deg(q_rel)
        
        difference_array[1:3, t] = euler_angles_deg
    end
    
    difference_array[4:6, :] = gt_target[5:7, :] - state_estimation_history[5:7, :]
    
    p = plot(layout=(6, 1), legend=false, size=(800, 1200))
    
    labels = ["Roll Difference", "Pitch Difference", "Yaw Difference", "Gyro Bias X", "Gyro Bias Y", "Gyro Bias Z"]
    
    for i in 1:6
        plot!(p[i], difference_array[i, :], title=labels[i], xlabel="Time", ylabel="Degrees")
    end
    
    display(p)
end