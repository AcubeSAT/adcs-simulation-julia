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
const norbits = 2
const qtarget = one(QuaternionF64)

const vecs = ADCSSims.generate_orbit_data(jd, norbits, 0.1)

const egm2008 = ADCSSims.GravityModels.load(ADCSSims.SatelliteToolboxGravityModels.IcgemFile, ADCSSims.SatelliteToolboxGravityModels.fetch_icgem_file(:EGM2008))
const n_max_dP = 1
P = Matrix{Float64}(undef, n_max_dP + 1, n_max_dP + 1)
dP = Matrix{Float64}(undef, n_max_dP + 1, n_max_dP + 1)

const PD = PDController(1e-2, 1e-1) # SMatrix{3,3}(I(3))
const qeci2body = ADCSSims.normalize(QuaternionF64(1))
const w = ADCSSims.@MVector [0.53, 0.53, 0.053]
sun_tracking = true
const state, τw, τsm, τgravs, τrmds = ADCSSims.rotational_dynamics(qeci2body, w, sun_tracking, PD, vecs..., qtarget, egm2008, n_max_dP, P, dP)
sun_tracking = false
const state2, τw2, τsm2, τgravs2, τrmds2 = ADCSSims.rotational_dynamics(state[end][2], state[end][1], sun_tracking, PD, vecs..., qtarget, egm2008, n_max_dP, P, dP)
sun_tracking = true
const state3, τw3, τsm3, τgravs3, τrmds3 = ADCSSims.rotational_dynamics(state2[end][2], state2[end][1], sun_tracking, PD, vecs..., qtarget, egm2008, n_max_dP, P, dP)


state = [state;state2;state3]

ADCSSims.plotτ(τw, τsm)
ADCSSims.plotwq(state)

q = [s[2] for s in state]
qorbit2body = [q1 * conj(q2) for (q1, q2) in zip(q, [vecs[7];vecs[7]])]
sun_eci = repeat(vecs[5],3)
nadir_eci = repeat(-normalize.(vecs[3]),3)

qbody2sun = [align_frame_with_vector(rotvec(sun_eci[i],q[i]), rotvec(nadir_eci[i],q[i]), [0,0,-1],[0,1,0]) for i in 1:length(q)]

ADCSSims.plotqs(qbody2sun)

ADCSSims.plotτgrav(τgravs)
ADCSSims.plotτgrav(τrmds)
