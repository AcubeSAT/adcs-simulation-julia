abstract type PointingMode end
struct SunPointing <: PointingMode
    q_offset::Quaternion
end
struct NadirPointing <: PointingMode
    q_offset::Quaternion
end
@concrete struct GroundTargetPointing <: PointingMode
    target_lat::Any
    target_lon::Any
    q_offset::Quaternion
end
@concrete struct PointingArguments
    sun_body::Any
    nadir_body::Any
    qeci2body::Any
    qeci2orbit::Any
    r_eci::Any
end

function mode_quaternion(Mode::PointingMode, args::PointingArguments)
    return error("Mode not implemented")
end
function mode_quaternion(Mode::SunPointing, args::PointingArguments)
    return Mode.q_offset *
           align_frame_with_vector(args.sun_body, args.nadir_body, [0, 0, -1], [0, 1, 0])
end
function mode_quaternion(Mode::NadirPointing, args::PointingArguments)
    return Mode.q_offset * (args.qeci2body * conj(args.qeci2orbit))
end
function mode_quaternion(Mode::GroundTargetPointing, args::PointingArguments)
    lat = deg2rad(Mode.target_lat)
    lon = deg2rad(Mode.target_lon)

    a = 6378.137  # Equatorial radius in kilometers
    b = 6356.752  # Polar radius in kilometers

    # Compute the radius at the given latitude using the WGS-84 model
    R_location = (a^2) / sqrt(a^2 * cos(lat)^2 + b^2 * sin(lat)^2)

    # Calculate the vector from the center of the Earth to the location using the computed radius
    x_location = R_location * cos(lat) * cos(lon)
    y_location = R_location * cos(lat) * sin(lon)
    z_location = R_location * sin(lat)
    vec_earth_to_location = [x_location, y_location, z_location]

    # Calculate the vector from the CubeSat to the location
    vec_cubesat_to_location_eci = normalize(vec_earth_to_location - args.r_eci)
    vec_cubesat_to_location_body = rotvec(vec_cubesat_to_location_eci, args.qeci2body)
    # Compute the quaternion for the desired rotation
    return Mode.q_offset * align_frame_with_vector(
        vec_cubesat_to_location_body, normalize(args.sun_body), [1.0, 0, 0], [0, 0, -1.0]
    )
end

function parse_pointing_mode(row)
    latitude = get(row, :latitude, missing)
    longitude = get(row, :longitude, missing)
    if row.pointing_mode == "NadirPointing"
        return NadirPointing(
            euler_to_quaternion(
                deg2rad(row.x_offset), deg2rad(row.y_offset), deg2rad(row.z_offset)
            ),
        )
    elseif row.pointing_mode == "GroundTargetPointing"
        return GroundTargetPointing(
            latitude,
            longitude,
            euler_to_quaternion(
                deg2rad(row.x_offset), deg2rad(row.y_offset), deg2rad(row.z_offset)
            ),
        )
    elseif row.pointing_mode == "SunPointing"
        return SunPointing(
            euler_to_quaternion(
                deg2rad(row.x_offset), deg2rad(row.y_offset), deg2rad(row.z_offset)
            ),
        )
    end
end
