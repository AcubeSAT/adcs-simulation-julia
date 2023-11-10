function calculate_orbit(jd, total_time, dt, orbital_elements)
    epc0 = Epoch(jd_to_caldate(jd)...)
    # oe0 = [R_EARTH + 522863.7, 0.01, 98.0, 306.615, 314.19, 99.89]
    eci0 = sOSCtoCART(orbital_elements, use_degrees = true)
    # T = orbit_period(oe0[1])
    epcf = epc0 + total_time
    orb = EarthInertialState(epc0, eci0, dt = dt, mass = 11.41, relativity = false)

    t, epc, eci = sim!(orb, epcf)
    return t, epc, eci
end

function rotational_dynamics(
    pointing_mode::PointingMode,
    t,
    epc::Vector{Epoch},
    r_eci,
    r_ecef,
    sun_eci,
    mag_eci,
    Qeci2orbit,
    R_ecef_to_eci,
    SimParams::SimulationParams,
    SimContext::SimulationContext,
    curindex,
)
    iters = length(t)
    println("iters from rotational_dynamics: $iters")
    t = epoch_to_datetime(epc)
    for i = 1:iters
        qeci2body = SimContext.state[curindex][2]
        nadir_body = -rotvec(normalize(r_eci[i]), qeci2body)
        sun_body = Vector(rotvec(sun_eci[i], qeci2body))
        mag_body = rotvec(mag_eci[i], qeci2body)
        target_vectors = (nadir_body, sun_body, sun_body)

        qeci2orbit = Qeci2orbit[i]
        PointingArgs =
            PointingArguments(sun_body, nadir_body, qeci2body, qeci2orbit, r_eci[i])

        wqeci2body, rτw, rτsm, τgrav, τrmd = control_loop(
            pointing_mode,
            SimParams,
            SimContext,
            PointingArgs,
            r_ecef[i],
            t[i],
            R_ecef_to_eci[i],
            mag_body,
            target_vectors,
            curindex,
        )

        curindex += 1
        SimContext.state[curindex] = wqeci2body
        SimContext.τw[curindex] = rτw
        SimContext.τsm[curindex] = rτsm
        SimContext.τgravs[curindex] = τgrav
        SimContext.τrmds[curindex] = τrmd
    end
    return curindex
end

function generate_orbit_data(jd, total_time, dt, orbital_elements)
    t, epc, eci = calculate_orbit(jd, total_time, dt, orbital_elements)
    r_eci = eci[1:3, :]
    r_eci = [vec(col) for col in eachcol(r_eci)]
    v_eci = eci[4:6, :]
    v_eci = [vec(col) for col in eachcol(v_eci)]
    rotation_eci2ecef = rECItoECEF.(epc)
    r_ecef = [rotation_eci2ecef[i] * r_eci[i] for i = 1:size(rotation_eci2ecef, 1)]
    rotation_ecef2eci = rECEFtoECI.(epc)
    year = jd_to_caldate(jd)[1]
    r, ϕ, θ = map(
        field -> [getfield(item, field) for item in SphericalFromCartesian().(r_ecef)],
        (:r, :ϕ, :θ),
    )
    mag_ecef = 1e-9 * ADCSSims.igrf.(year, r, ϕ, θ)
    mag_eci = [rotation_ecef2eci[i] * mag_ecef[i] for i = 1:size(rotation_ecef2eci, 1)]
    sun_eci = sun_position.(epc)
    sun_eci = normalize.(sun_eci)
    T = eci2orbit.(r_eci, v_eci)
    qeci2orbit = from_rotation_matrix.(T)
    return t, epc, r_eci, r_ecef, sun_eci, mag_eci, qeci2orbit, rotation_ecef2eci
end

function eci2orbit(r_eci, v_eci)
    z = -normalize(r_eci)
    y = normalize(cross(r_eci, v_eci))
    x = cross(y, z)
    T = hcat(x, y, z)
    return T'
end
