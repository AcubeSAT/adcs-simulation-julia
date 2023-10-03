using ConcreteStructs, LinearAlgebra, StaticArrays, Plots

function plotwq(state)
    q1 = [s[2].coeffs[1] for s in state]
    q2 = [s[2].coeffs[2] for s in state]
    q3 = [s[2].coeffs[3] for s in state]
    q4 = [s[2].coeffs[4] for s in state]

    v1 = [s[1][1] for s in state]
    v2 = [s[1][2] for s in state]
    v3 = [s[1][3] for s in state]

    p1 = plot(q1, label="q1", title="Quaternion q1")
    p2 = plot(q2, label="q2", title="Quaternion q2")
    p3 = plot(q3, label="q3", title="Quaternion q3")
    p4 = plot(q4, label="q4", title="Quaternion q4")
    p5 = plot(v1, label="v1", title="Vector v1")
    p6 = plot(v2, label="v2", title="Vector v2")
    p7 = plot(v3, label="v3", title="Vector v3")

    plt = plot(p1, p2, p3, p4, p5, p6, p7, layout=(3,3), legend=false)

    display(plt)
    return nothing
end

function plotτs(state)
    p5 = plot(state[1, :], label="v1", title="Vector v1")
    p6 = plot(state[2, :], label="v2", title="Vector v2")
    p7 = plot(state[3, :], label="v3", title="Vector v3")

    plt = plot(p5, p6, p7, layout=(1,3), legend=false)

    display(plt)
    return nothing
end

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

@concrete struct PDController
    Kp
    Kd
end

function calculate_torque(PD::PDController, qtarget, qestimated, w, wtarget)
    qrel = qestimated * conj(qtarget)
    werr = w - wtarget
    return -sign(real(qrel)) * PD.Kp * vec(qrel) - PD.Kd * werr
end

function control(iters, dt)
    q = Quaternion(0.1531, 0.6853, 0.6953, 0.1531)
    qtarget = normalize(QuaternionF64(1,-2,3,6))
    qestimated = normalize(Quaternion(q.coeffs .+ 1e-3))
    w = @MVector [0.53, 0.53, 0.053]
    wtarget = zeros(3)
    Inertia = diagm([10000,9000,12000])
    τs = Vector{MVector{3, Float64}}(undef, iters)
    res = Vector{Tuple{typeof(w), typeof(q)}}(undef, iters)
    PD = PDController(50, 500)
    for i in 1:iters
        τ = calculate_torque(PD, qtarget, qestimated, w, wtarget)
        τs[i] = τ
        w, q = rk4(Inertia, w, τ, q, dt)
        qestimated = normalize(Quaternion(q.coeffs .+ 1e-3))
        res[i] = w, q
    end
    return res, τs
end

res, τs = control(3000, 0.1);
plotwq(res)
plotτs(reduce(hcat, τs))
