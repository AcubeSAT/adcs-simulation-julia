using ADCSSims
using BenchmarkTools
using JET
using DataFrames
using CSV
using TOML

function get_total_time(df::DataFrame)
    sum(skipmissing(df[!, "duration"]))
end

const config = TOML.parsefile("config.toml")
const schedule_df = CSV.File("schedule.csv", types=[Float64, String, Union{Missing, Float64}, Union{Missing, Float64}]) |> DataFrame

const jd = config["simulation"]["jd"]
const qtarget = ADCSSims.Quaternion(config["simulation"]["qtarget"])
const dt = config["simulation"]["dt"]
const total_time = get_total_time(schedule_df)
const vecs = ADCSSims.generate_orbit_data(jd, total_time, dt)
const egm2008 = ADCSSims.GravityModels.load(ADCSSims.SatelliteToolboxGravityModels.IcgemFile, ADCSSims.SatelliteToolboxGravityModels.fetch_icgem_file(:EGM2008))
const n_max_dP = 1
const P = Matrix{Float64}(undef, n_max_dP + 1, n_max_dP + 1)
const dP = Matrix{Float64}(undef, n_max_dP + 1, n_max_dP + 1)
const PD = PDController(0.2, 1.0)

function run_pointing_modes(df::DataFrame)
    cumulative_start_time = 0.0
    qeci2body = one(QuaternionF64)
    w = ADCSSims.@MVector [0.0, 0.0, 0.0]
    state_history = []  
    τw_history, τsm_history, τgravs_history, τrmds_history = [], [], [], []
    for row in eachrow(df)
        start_time = cumulative_start_time
        end_time = start_time + row.duration
        pointing_mode = ADCSSims.parse_pointing_mode(row.pointing_mode)
        latitude = get(row, :latitude, missing)
        longitude = get(row, :longitude, missing) 
        start_index = Int(floor(start_time / dt)) + 1
        end_index = Int(ceil(end_time / dt))
        vectors_slice = ADCSSims.subvector(vecs, start_index, end_index)
        state, τw, τsm, τgravs, τrmds = ADCSSims.rotational_dynamics(
            qeci2body, w, pointing_mode, ADCSSims.StaticPointingArgs(latitude, longitude), PD, vectors_slice..., dt, qtarget, egm2008, n_max_dP, P, dP
        )
        append!(state_history, state)
        append!(τw_history, τw)
        append!(τsm_history, τsm)
        append!(τgravs_history, τgravs)
        append!(τrmds_history, τrmds)

        qeci2body = state[end][2]  
        w = state[end][1]  
        println("From $start_time to $end_time, mode: $pointing_mode")
        cumulative_start_time = end_time
    end
    return state_history, τw_history, τsm_history, τgravs_history, τrmds_history
end

state_history, τw_history, τsm_history, τgravs_history, τrmds_history = run_pointing_modes(schedule_df)
# 5705.307041952439 total period
# 21                capture image (nadir for now)
# 479.6297          GS tracking (nadir for now)
# 5204.677341952439 sun tracking

q = [s[2] for s in state_history]
qorbit2body = [q1 * conj(q2) for (q1, q2) in zip(q, [vecs[7]; vecs[7]])]
sun_eci = vecs[5]
nadir_eci = -ADCSSims.normalize.(vecs[3])

qbody2sun = [align_frame_with_vector(rotvec(sun_eci[i], q[i]), rotvec(nadir_eci[i], q[i]), [0, 0, -1], [0, 1, 0]) for i in 1:length(q)]

ADCSSims.plotqs(qbody2sun)
# ADCSSims.plotτgrav(τgravs)
# ADCSSims.plotτgrav(τrmds)
# ADCSSims.plotτ(τw_history, τsm_history)
# ADCSSims.plotwq(state_history)

n = 1

jd_values = [ADCSSims.jd(epc) for epc in vecs[2][1:n:end]]

coeff1 = [q[1] for q in qbody2sun[1:n:end]]
coeff2 = [q[2] for q in qbody2sun[1:n:end]]
coeff3 = [q[3] for q in qbody2sun[1:n:end]]
coeff4 = [q[4] for q in qbody2sun[1:n:end]]

DataFrame(
    JD=jd_values,
    q1=coeff1,
    q2=coeff2,
    q3=coeff3,
    q4=coeff4
) |> CSV.write("data.csv")
