@concrete struct PDController
    Kp
    Kd
end

function calculate_torque(PD::PDController, qtarget, qestimated, w, wtarget, qeci2body)
    qrel = qestimated * conj(qtarget)
    w_io_i = @SVector [0.0, -0.00014277091387915417, 0.0010992751451387854] # w of orbit frame, rad/s
    w_io_b = rotvec(w_io_i, qeci2body)
    wrel = w - w_io_b
    werr = wrel - wtarget
    return -sign(real(qrel)) * PD.Kp * vec(qrel) - PD.Kd * werr
end

function decompose_torque(τ, b, msaturation)
    τw = dot(τ, b) / dot(b, b) * b
    τm = τ - τw
    m = cross(b, τm) / dot(b, b)
    mprime = m / norm(τm)
    ksm = min(
        abs(msaturation / mprime[1]),
        abs(msaturation / mprime[2]),
        abs(msaturation / mprime[3]),
    )
    mtrue = ksm * mprime
    τsm = ksm * cross(b, mprime)
    τw = τ - τsm
    return τw, τsm, mtrue
end

# TODO: what if saturation compensation is smaller than the cubesat w from control
# TODO: reaction wheel not taken into account for anything
function control_loop(
    Mode::PointingMode,
    SimParams::SimulationParams,
    SimContext::SimulationContext,
    PointingArgs::PointingArguments,
    mag_body,
    target_vectors,
    curindex,
)
    w, qeci2body = SimContext.state[curindex]
    qerr = emulate_estimation(SimParams.sensors, target_vectors, w)
    qestimated = qerr * mode_quaternion(Mode, PointingArgs)
    τ = calculate_torque(
        SimParams.PD, SimParams.qtarget, qestimated, w, SimParams.wtarget, qeci2body
    )
    τw, τsm, mtrue = decompose_torque(τ, mag_body, SimParams.msaturation)
    compensation =
        deadzone_compensation(SimContext.RW) + saturation_compensation(SimContext.RW)
    τw = clamp.(τw + compensation, -SimContext.RW.maxtorque, SimContext.RW.maxtorque)
    # rwfriction = stribeck(RW)
    # @reset SimContext.RW.w = rk4_rw(SimContext.RW.J, SimContext.RW.w, -τw, SimParams.dt)
    return τw, τsm
end
