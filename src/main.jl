using ADCSSims

using Optim

inertia_matrix = ADCSSims.diagm([1.0; 1.5; 0.5])
sigma_u = 7.7570e-5;
sigma_v = 0.0026;
mag_noise = 1e-3
sun_noise = 1e-3
dt = 0.1
constant_params = ADCSSims.parameter_struct(inertia_matrix,
    sigma_u,
    sigma_v,
    mag_noise,
    sun_noise,
    dt)

(q_history, w_history, bias_history, mag_eci, sun_eci, mag_noisy_history, sun_noisy_history, gyro_noisy_history) = ADCSSims.run_groundtruth_simulation(constant_params)
sun_noisy_history[:, 5001:6500] = zeros(3, 1500)

function objective_function(x)
    tunable_params = ADCSSims.package_weights(x)
    statepred = ADCSSims.run_filter_simulation(tunable_params,
        constant_params,
        mag_noisy_history,
        sun_noisy_history,
        mag_eci,
        sun_eci,
        gyro_noisy_history)

    q̂ = [sp.q for sp in statepred]
    return sum(ADCSSims.qloss.(q̂, q_history))
end

initial_x = [-6.0, -16, -6, -2]
result = optimize(objective_function,
    initial_x,
    Optim.NelderMead(),
    Optim.Options(iterations = 200, show_trace = true, g_tol = 1e-15))

optimized_x = Optim.minimizer(result)
tunable_params = ADCSSims.package_weights(optimized_x)
ŷ = ADCSSims.run_filter_simulation(tunable_params,
    constant_params,
    mag_noisy_history,
    sun_noisy_history,
    mag_eci,
    sun_eci,
    gyro_noisy_history)

len = length(ŷ)
y = Vector{ADCSSims.KFState}(undef, len)
for i in 1:len
    y[i] = ADCSSims.KFState(q_history[i], bias_history[:, i])
end
ADCSSims.plot_histories(ŷ, y)
ADCSSims.plot_difference(y, ŷ)

const jd = 2459921.0
const norbits = 4
const qtarget = Quaternion(1.0,0.0,0.0,0.0)
const vecs = ADCSSims.generate_orbit_data(jd, norbits, 0.1)
const PD = PDController(1e-4, 2e-4) # SMatrix{3,3}(I(3))
const state, τw, τsm = ADCSSims.rotational_dynamics(PD, vecs..., qtarget)

plotτ(τw, τsm)
plotqs(state)

function plotτ(τw, τsm)
    τw1 = [x[1] for x in τw]
    τw2 = [x[2] for x in τw]
    τw3 = [x[3] for x in τw]
    τsm1 = [x[1] for x in τsm]
    τsm2 = [x[2] for x in τsm]
    τsm3 = [x[3] for x in τsm]

    p1 = plot(τw1, label="τw1", title="RW torque 1")
    p2 = plot(τw2, label="τw2", title="RW torque 2")
    p3 = plot(τw3, label="τw3", title="RW torque 3")
    p4 = plot(τsm1, label="τsm1", title="MTQ torque 1")
    p5 = plot(τsm2, label="τsm2", title="MTQ torque 2")
    p6 = plot(τsm3, label="τsm3", title="MTQ torque 3")

    plt = plot(p1, p2, p3, p4, p5, p6, layout=(2,3), legend=false, size=(700,300))
    display(plt)
    return nothing
end

function plotqs(state)
    # Extracting quaternion components and vector elements
    q1 = [s[2].coeffs[1] for s in state]
    q2 = [s[2].coeffs[2] for s in state]
    q3 = [s[2].coeffs[3] for s in state]
    q4 = [s[2].coeffs[4] for s in state]

    v1 = [s[1][1] for s in state]
    v2 = [s[1][2] for s in state]
    v3 = [s[1][3] for s in state]

    # Creating subplots
    p1 = plot(q1, label="q1", title="Quaternion q1")
    p2 = plot(q2, label="q2", title="Quaternion q2")
    p3 = plot(q3, label="q3", title="Quaternion q3")
    p4 = plot(q4, label="q4", title="Quaternion q4")
    p5 = plot(v1, label="v1", title="Vector v1")
    p6 = plot(v2, label="v2", title="Vector v2")
    p7 = plot(v3, label="v3", title="Vector v3")

    # Combining subplots into a single plot with a 3x3 grid
    # The 8th and 9th subplots remain empty
    plt = plot(p1, p2, p3, p4, p5, p6, p7, layout=(3,3), legend=false)

    # Display the plots (if you're in a REPL environment or Jupyter Notebook)
    display(plt)
    return nothing
end

