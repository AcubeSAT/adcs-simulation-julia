using ADCSSims
using DataFrames
using CSV
using Base.Iterators: partition

const jd = 2459921.0
const norbits = 1
const qtarget = one(QuaternionF64)
const dt = 0.1

const vecs = ADCSSims.generate_orbit_data(jd, norbits, dt)

function split_into_parts(vec::Vector, n::Int)
    len = length(vec)
    size_per_part = div(len, n)

    parts = Vector{typeof(vec)}(undef, n)
    start_idx = 1

    for i in 1:n-1
        parts[i] = vec[start_idx:start_idx+size_per_part-1]
        start_idx += size_per_part
    end
    parts[n] = vec[start_idx:end] # last part takes the remainder

    return parts
end

const split_vecs = [split_into_parts(vec, 3) for vec in vecs]

const first_parts = [v[1] for v in split_vecs]
const second_parts = [v[2] for v in split_vecs]
const third_parts = [v[3] for v in split_vecs]

const egm2008 = ADCSSims.GravityModels.load(ADCSSims.SatelliteToolboxGravityModels.IcgemFile, ADCSSims.SatelliteToolboxGravityModels.fetch_icgem_file(:EGM2008))
const n_max_dP = 1
P = Matrix{Float64}(undef, n_max_dP + 1, n_max_dP + 1)
dP = Matrix{Float64}(undef, n_max_dP + 1, n_max_dP + 1)

const PD = PDController(1e-2, 1e-1) # SMatrix{3,3}(I(3))
const qeci2body = one(QuaternionF64)
const w = ADCSSims.@MVector [0.0, 0.0, 0.0]
const state, τw, τsm, τgravs, τrmds = ADCSSims.rotational_dynamics(qeci2body, w, SunPointing, PD, first_parts..., dt, qtarget, egm2008, n_max_dP, P, dP)
const state2, τw2, τsm2, τgravs2, τrmds2 = ADCSSims.rotational_dynamics(state[end][2], state[end][1], NadirPointing, PD, second_parts..., dt, qtarget, egm2008, n_max_dP, P, dP)
const state3, τw3, τsm3, τgravs3, τrmds3 = ADCSSims.rotational_dynamics(state2[end][2], state2[end][1], SunPointing, PD, third_parts..., dt, qtarget, egm2008, n_max_dP, P, dP)

state_full = [state; state2; state3]

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
