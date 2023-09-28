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
const norbits = 1
const qtarget = one(QuaternionF64)
const vecs = ADCSSims.generate_orbit_data(jd, norbits, 0.1)
const PD = PDController(1e-4, 1e-3) # SMatrix{3,3}(I(3))
const state, τw, τsm = ADCSSims.rotational_dynamics(PD, vecs..., qtarget)

ADCSSims.plotτ(τw, τsm)
ADCSSims.plotqs(state)
