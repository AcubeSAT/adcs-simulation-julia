function plot_histories(ŷ, y)
    get_comp(comp_func, array) = [comp_func(s) for s in array]
    funcs = [
        s -> s.q.q1,
        s -> s.q.q2,
        s -> s.q.q3,
        s -> s.q.q4,
        s -> s.bias[1],
        s -> s.bias[2],
        s -> s.bias[3],
    ]
    titles = [
        "q1 values",
        "q2 values",
        "q3 values",
        "q4 values",
        "bias1 values",
        "bias2 values",
        "bias3 values",
    ]
    p = plot(layout = (7, 1), size = (600, 1200))
    for i in 1:7
        ŷ_vals = get_comp(funcs[i], ŷ)
        y_vals = get_comp(funcs[i], y)
        plot!(p[i],
            ŷ_vals,
            label = "Estimated",
            xlabel = "Timestep",
            ylabel = "Value",
            linewidth = 2,
            title = titles[i])
        plot!(p[i],
            y_vals,
            label = "Ground Truth",
            xlabel = "Timestep",
            ylabel = "Value",
            linewidth = 2)
    end
    display(p)
end

function plot_difference(ŷ, y)
    len = length(y)
    difference_array = zeros(Float64, 6, len)
    ŷqs = [s.q for s in ŷ]
    yqs = [s.q for s in y]
    for i in 1:len
        q_rel = ŷqs[i] * conj(yqs[i])
        euler_angles_deg = quat_to_euler_deg(q_rel)
        difference_array[1:3, i] = euler_angles_deg
    end

    ŷbias = reduce(hcat, [s.bias for s in ŷ])
    ybias = reduce(hcat, [s.bias for s in y])
    difference_array[4:6, :] = ybias - ŷbias
    p = plot(layout = (6, 1), legend = false, size = (800, 1200))
    labels = [
        "Roll Difference",
        "Pitch Difference",
        "Yaw Difference",
        "Gyro Bias X",
        "Gyro Bias Y",
        "Gyro Bias Z",
    ]

    for i in 1:6
        plot!(p[i],
            difference_array[i, :],
            title = labels[i],
            xlabel = "Time",
            ylabel = "Degrees")
    end
    display(p)
end
