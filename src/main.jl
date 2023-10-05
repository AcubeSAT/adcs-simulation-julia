using ADCSSims
using DataFrames
using CSV
using Base.Iterators: partition

const jd = 2459921.1
const norbits = 1
const qtarget = one(QuaternionF64)
const dt = 0.1

const vecs = ADCSSims.generate_orbit_data(jd, norbits, dt)
# 5705.307041952439 total period
# 21                capture image (nadir for now)
# 479,6297          GS tracking (nadir for now)
# 5204.677341952439 sun tracking
const vsunpre = ADCSSims.subvector(vecs, 1, 26020)
const vimage = ADCSSims.subvector(vecs, 26021, 26231)
const vgs = ADCSSims.subvector(vecs, 26232, 31028)
const vsunpost = ADCSSims.subvector(vecs, 31029, :end)

const egm2008 = ADCSSims.GravityModels.load(ADCSSims.SatelliteToolboxGravityModels.IcgemFile, ADCSSims.SatelliteToolboxGravityModels.fetch_icgem_file(:EGM2008))
const n_max_dP = 1
P = Matrix{Float64}(undef, n_max_dP + 1, n_max_dP + 1)
dP = Matrix{Float64}(undef, n_max_dP + 1, n_max_dP + 1)

const PD = PDController(1e-2, 1e-1) # SMatrix{3,3}(I(3))
const qeci2body = one(QuaternionF64)
const w = ADCSSims.@MVector [0.0, 0.0, 0.0]
const state, τw, τsm, τgravs, τrmds = ADCSSims.rotational_dynamics(qeci2body, w, SunPointing, PD, vsunpre..., dt, qtarget, egm2008, n_max_dP, P, dP)
const state2, τw2, τsm2, τgravs2, τrmds2 = ADCSSims.rotational_dynamics(state[end][2], state[end][1], NadirPointing, PD, vimage..., dt, qtarget, egm2008, n_max_dP, P, dP)
const state3, τw3, τsm3, τgravs3, τrmds3 = ADCSSims.rotational_dynamics(state2[end][2], state2[end][1], NadirPointing, PD, vgs..., dt, qtarget, egm2008, n_max_dP, P, dP)
const state4, τw4, τsm4, τgravs4, τrmds4 = ADCSSims.rotational_dynamics(state3[end][2], state3[end][1], SunPointing, PD, vsunpost..., dt, qtarget, egm2008, n_max_dP, P, dP)

state_full = [state; state2; state3; state4]

ADCSSims.plotτ(τw, τsm)
ADCSSims.plotwq(state_full)

q = [s[2] for s in state_full]
qorbit2body = [q1 * conj(q2) for (q1, q2) in zip(q, [vecs[7]; vecs[7]])]
sun_eci = vecs[5]
nadir_eci = -ADCSSims.normalize.(vecs[3])

qbody2sun = [align_frame_with_vector(rotvec(sun_eci[i], q[i]), rotvec(nadir_eci[i], q[i]), [0, 0, -1], [0, 1, 0]) for i in 1:length(q)]

ADCSSims.plotqs(qbody2sun)

ADCSSims.plotτgrav(τgravs)
ADCSSims.plotτgrav(τrmds)

n = 100

jd_values = [ADCSSims.jd(epc) for epc in vecs[2][1:n:end]]

coeff1 = [q.coeffs[1] for q in qbody2sun[1:n:end]]
coeff2 = [q.coeffs[2] for q in qbody2sun[1:n:end]]
coeff3 = [q.coeffs[3] for q in qbody2sun[1:n:end]]
coeff4 = [q.coeffs[4] for q in qbody2sun[1:n:end]]

DataFrame(
    JD=jd_values,
    q1=coeff1,
    q2=coeff2,
    q3=coeff3,
    q4=coeff4
) |> CSV.write("data.csv")
