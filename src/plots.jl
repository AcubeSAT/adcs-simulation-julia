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
    for i = 1:7
        ŷ_vals = get_comp(funcs[i], ŷ)
        y_vals = get_comp(funcs[i], y)
        plot!(
            p[i],
            ŷ_vals,
            label = "Estimated",
            xlabel = "Timestep",
            ylabel = "Value",
            linewidth = 2,
            title = titles[i],
        )
        plot!(
            p[i],
            y_vals,
            label = "Ground Truth",
            xlabel = "Timestep",
            ylabel = "Value",
            linewidth = 2,
        )
    end
    display(p)
    return nothing
end

function plot_difference(ŷ, y)
    len = length(y)
    difference_array = zeros(Float64, 6, len)
    ŷqs = [s.q for s in ŷ]
    yqs = [s.q for s in y]
    for i = 1:len
        q_rel = ŷqs[i] * conj(yqs[i])
        euler_angles_deg = toeuler(q_rel)
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

    for i = 1:6
        plot!(
            p[i],
            difference_array[i, :],
            title = labels[i],
            xlabel = "Time",
            ylabel = "Degrees",
        )
    end
    display(p)
    return nothing
end


function plotτ(τw, τsm)
    τw1 = [x[1] for x in τw]
    τw2 = [x[2] for x in τw]
    τw3 = [x[3] for x in τw]
    τsm1 = [x[1] for x in τsm]
    τsm2 = [x[2] for x in τsm]
    τsm3 = [x[3] for x in τsm]

    p1 = plot(τw1, label = "τw1", title = "RW torque 1")
    p2 = plot(τw2, label = "τw2", title = "RW torque 2")
    p3 = plot(τw3, label = "τw3", title = "RW torque 3")
    p4 = plot(τsm1, label = "τsm1", title = "MTQ torque 1")
    p5 = plot(τsm2, label = "τsm2", title = "MTQ torque 2")
    p6 = plot(τsm3, label = "τsm3", title = "MTQ torque 3")

    plt = plot(p1, p2, p3, p4, p5, p6, layout = (2, 3), legend = false, size = (700, 300))
    display(plt)
    return nothing
end

function plotτgrav(τgravs)
    τg1 = [x[1] for x in τgravs]
    τg2 = [x[2] for x in τgravs]
    τg3 = [x[3] for x in τgravs]

    p1 = plot(τg1, label = "τg1", title = "Gravity torque 1")
    p2 = plot(τg2, label = "τg2", title = "Gravity torque 2")
    p3 = plot(τg3, label = "τg3", title = "Gravity torque 3")

    plt = plot(p1, p2, p3, layout = (3, 1), legend = false, size = (700, 300))
    display(plt)
    return nothing
end

function plotwq(state)
    # Extracting quaternion components and vector elements
    q1 = [s[2].coeffs[1] for s in state]
    q2 = [s[2].coeffs[2] for s in state]
    q3 = [s[2].coeffs[3] for s in state]
    q4 = [s[2].coeffs[4] for s in state]

    v1 = [s[1][1] for s in state]
    v2 = [s[1][2] for s in state]
    v3 = [s[1][3] for s in state]

    p1 = plot(q1, label = "q1", title = "Quaternion q1")
    p2 = plot(q2, label = "q2", title = "Quaternion q2")
    p3 = plot(q3, label = "q3", title = "Quaternion q3")
    p4 = plot(q4, label = "q4", title = "Quaternion q4")
    p5 = plot(v1, label = "v1", title = "Vector v1")
    p6 = plot(v2, label = "v2", title = "Vector v2")
    p7 = plot(v3, label = "v3", title = "Vector v3")

    plt = plot(p1, p2, p3, p4, p5, p6, p7, layout = (3, 3), legend = false)

    display(plt)
    return nothing
end

function plotqs(qorbit2body)
    q1 = [q.coeffs[1] for q in qorbit2body]
    q2 = [q.coeffs[2] for q in qorbit2body]
    q3 = [q.coeffs[3] for q in qorbit2body]
    q4 = [q.coeffs[4] for q in qorbit2body]

    p1 = plot(q1, label = "q1", title = "Quaternion q1")
    p2 = plot(q2, label = "q2", title = "Quaternion q2")
    p3 = plot(q3, label = "q3", title = "Quaternion q3")
    p4 = plot(q4, label = "q4", title = "Quaternion q4")

    plt = plot(p1, p2, p3, p4, layout = (2, 2), legend = false)

    display(plt)
    return nothing
end
