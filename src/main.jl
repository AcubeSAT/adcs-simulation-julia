using ADCSSims

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
const state, τw, τsm, τgravs, τrmds = ADCSSims.rotational_dynamics(qeci2body, w, SunPointing, PD, vecs..., qtarget, egm2008, n_max_dP, P, dP)
const state2, τw2, τsm2, τgravs2, τrmds2 = ADCSSims.rotational_dynamics(state[end][2], state[end][1], NadirPointing, PD, vecs..., qtarget, egm2008, n_max_dP, P, dP)
const state3, τw3, τsm3, τgravs3, τrmds3 = ADCSSims.rotational_dynamics(state2[end][2], state2[end][1], SunPointing, PD, vecs..., qtarget, egm2008, n_max_dP, P, dP)

state_full = [state;state2;state3]

ADCSSims.plotτ(τw, τsm)
ADCSSims.plotwq(state_full)

q = [s[2] for s in state_full]
qorbit2body = [q1 * conj(q2) for (q1, q2) in zip(q, [vecs[7];vecs[7]])]
sun_eci = repeat(vecs[5],3)
nadir_eci = repeat(-ADCSSims.normalize.(vecs[3]),3)

qbody2sun = [align_frame_with_vector(rotvec(sun_eci[i],q[i]), rotvec(nadir_eci[i],q[i]), [0,0,-1],[0,1,0]) for i in 1:length(q)]

ADCSSims.plotqs(qbody2sun)

ADCSSims.plotτgrav(τgravs)
ADCSSims.plotτgrav(τrmds)
