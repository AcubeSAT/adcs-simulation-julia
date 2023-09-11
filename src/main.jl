using SatelliteToolboxGeomagneticField
using SatelliteDynamics
using Flux
using Optim

include("dynamics.jl")
include("utils.jl")
include("parameters.jl")
include("mekf.jl")
include("noise_models.jl")
include("simulation.jl")
include("plots.jl")

inertia_matrix = diagm([1.0;1.5;0.5])
sigma_u = 7.7570e-5;
sigma_v = 0.0026;
mag_noise = 1e-3
sun_noise = 1e-4
dt = 0.1
constant_params = parameter_struct(inertia_matrix, sigma_u, sigma_v, mag_noise, sun_noise, dt)

(q_history, w_history, bias_history, mag_eci, sun_eci, mag_noisy_history, sun_noisy_history, gyro_noisy_history) = run_groundtruth_simulation(constant_params)
gt_target = [q_history;bias_history]

# Make eclipse
sun_noisy_history[:, 5001:6500] = zeros(3,1500)

function objective_function(x)
    tunable_params = package_weights(x)
    state_estimation_history = run_filter_simulation(tunable_params, constant_params, mag_noisy_history, sun_noisy_history, mag_eci, sun_eci, gyro_noisy_history)
    loss = Flux.mse(gt_target[1:4,:], state_estimation_history[1:4,:])
end

# # Initial guess
initial_x = [-9,-10.6,-6.5,-13]

# Perform optimization using Nelder-Mead method
result = optimize(objective_function, initial_x, NelderMead(), Optim.Options(iterations = 100, show_trace=true, g_tol=1e-15))

# Extract optimized parameters
optimized_x = Optim.minimizer(result)

tunable_params = package_weights(optimized_x)
state_estimation_history = run_filter_simulation(tunable_params, constant_params, mag_noisy_history, sun_noisy_history, mag_eci, sun_eci, gyro_noisy_history)
plot_histories(gt_target, state_estimation_history) 
plot_difference(gt_target, state_estimation_history)