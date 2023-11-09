using ADCSSims
using BenchmarkTools
using JET
using DataFrames
using CSV
using TOML

function get_total_time(df::DataFrame)
    sum(skipmissing(df[!, "duration"]))
end

ADCSSims.@concrete struct ConfigParams
    jd
    dt
    qtarget
    wtarget
    w
    Kp
    Kd
    grmodel
    n_max_dP
    semi_major_axis
    eccentricity
    inclination
    RAAN
    argument_of_perigee
    mean_anomaly
    rmd
    msaturation
    total_time
end

# orbital_elements = [ADCSSims.R_EARTH + 522863.7, 0.01, 98.0, 306.615, 314.19, 99.89]
function parseconfig()
    config = TOML.parsefile("config.toml")
    schedule_df = CSV.File("schedule.csv", types=[Float64, String, Union{Missing,Float64}, Union{Missing,Float64}]) |> DataFrame

    return schedule_df, ConfigParams(config["simulation"]["jd"],
        config["simulation"]["dt"],
        ADCSSims.Quaternion(config["simulation"]["qtarget"]),
        config["simulation"]["wtarget"],
        config["simulation"]["w"],
        config["controller"]["Kp"],
        config["controller"]["Kd"],
        Symbol(config["gravity"]["model"]),
        config["gravity"]["n_max_dP"],
        config["orbit"]["semi_major_axis"] + ADCSSims.R_EARTH,
        config["orbit"]["eccentricity"],
        config["orbit"]["inclination"],
        config["orbit"]["RAAN"],
        config["orbit"]["argument_of_perigee"],
        config["orbit"]["mean_anomaly"],
        ADCSSims.SVector{3}(config["disturbances"]["rmd"]),
        config["actuators"]["msaturation"],
        get_total_time(schedule_df))
end

function init()
    schedule_df, config = parseconfig()
    vecs = ADCSSims.generate_orbit_data(config.jd, config.total_time, config.dt,
        [config.semi_major_axis,
            config.eccentricity,
            config.inclination,
            config.RAAN,
            config.argument_of_perigee,
            config.mean_anomaly])

    egm2008 = ADCSSims.GravityModels.load(
        ADCSSims.SatelliteToolboxGravityModels.IcgemFile,
        ADCSSims.SatelliteToolboxGravityModels.fetch_icgem_file(config.grmodel))

    P = Matrix{Float64}(undef, config.n_max_dP + 1, config.n_max_dP + 1)
    dP = Matrix{Float64}(undef, config.n_max_dP + 1, config.n_max_dP + 1)

    PD = PDController(config.Kp, config.Kd)
    sensors = (NadirSensor(), StarTracker(), SunSensor())

    Ixx = 0.228128
    Iyy = 0.248027
    Izz = 0.091558
    Ixy = 0.000147
    Iyz = -0.000086
    Izx = 0.024513
    inertia_matrix = ADCSSims.@SMatrix [[Ixx, Ixy, Izx] [Ixy, Iyy, Iyz] [Izx, Iyz, Izz]]
    SimParams = SimulationParams(PD,
        config.qtarget,
        config.wtarget,
        egm2008,
        config.n_max_dP,
        P,
        dP,
        config.msaturation,
        sensors,
        inertia_matrix,
        config.dt,
        config.rmd)

    niter = length(vecs[1]) + 1
    state_history = Vector{Tuple{Vector{Float64},QuaternionF64}}(undef, niter)
    state_history[1] = (config.w, one(QuaternionF64))
    τw = Vector{Vector{Float64}}(undef, niter)
    τsm = Vector{Vector{Float64}}(undef, niter)
    τgrav = Vector{Vector{Float64}}(undef, niter)
    τrmd = Vector{Vector{Float64}}(undef, niter)

    RW = ReactionWheel(J=ADCSSims.I(3),
        w=94.247779 * ones(3),
        saturationα=1,
        deadzoneα=1,
        maxtorque=0.001)

    SimContext = SimulationContext(state_history,
        τw,
        τsm,
        τgrav,
        τrmd,
        vecs[6],
        RW)

    curindex = 1

    return SimParams, SimContext, schedule_df, vecs, curindex
end

function run_pointing_modes(SimParams::SimulationParams, SimContext::SimulationContext, df::DataFrame, vecs, curindex)
    cumulative_start_time = 0.0
    for row in eachrow(df)
        start_time = cumulative_start_time
        end_time = start_time + row.duration
        pointing_mode = ADCSSims.parse_pointing_mode(row)
        start_index = Int(floor(start_time / SimParams.dt)) + 1
        end_index = Int(ceil(end_time / SimParams.dt))
        vectors_slice = ADCSSims.subvector(vecs, start_index, end_index)
        curindex = ADCSSims.rotational_dynamics(pointing_mode, vectors_slice..., SimParams, SimContext, curindex)
        println("From $start_time to $end_time, mode: $pointing_mode")
        cumulative_start_time = end_time
    end
    return SimContext
end

function main()
    SimParams, SimContext, schedule_df, vecs, curindex = init()
    run_pointing_modes(SimParams, SimContext, schedule_df, vecs, curindex)

    println("sc length: $(lastindex(SimContext.state))")
    q = [s[2] for s in SimContext.state]
    ADCSSims.plotqs(q)
    sun_eci = vecs[5]
    nadir_eci = -ADCSSims.normalize.(vecs[3])

    qbody2sun = [align_frame_with_vector(rotvec(sun_eci[i], q[i]), rotvec(nadir_eci[i], q[i]), [0, 0, -1], [0, 1, 0]) for i in 1:length(q)-1]
    ADCSSims.plotqs(qbody2sun)

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
end

# 5705.307041952439 total period
# 21                capture image (nadir for now)
# 479.6297          GS tracking (nadir for now)
# 5204.677341952439 sun tracking
