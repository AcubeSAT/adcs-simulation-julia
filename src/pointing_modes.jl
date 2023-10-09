abstract type PointingMode end
struct SunPointing <: PointingMode end
struct NadirPointing <: PointingMode end
struct GSPointing <: PointingMode end
struct CamPointing <: PointingMode end

@concrete struct StaticPointingArgs
    gs_lat
    gs_lon
    cam_lat
    cam_lon
end

@concrete struct DynamicPointingArgs
    sun_body
    nadir_body
    qeci2body
    qeci2orbit
    r_eci
end

@concrete struct PointingArguments
    dynamic_args::DynamicPointingArgs
    static_args::StaticPointingArgs
end

function mode_quaternion(::Type{PointingMode}, args::PointingArguments)
    error("Mode not implemented")
end
function mode_quaternion(::Type{SunPointing}, args::PointingArguments)
    return align_frame_with_vector(args.dynamic_args.sun_body, args.dynamic_args.nadir_body, [0,0,-1], [0,1,0])
end
function mode_quaternion(::Type{NadirPointing}, args::PointingArguments)
    return args.dynamic_args.qeci2body * conj(args.dynamic_args.qeci2orbit)
end
function mode_quaternion(::Type{GSPointing}, args::PointingArguments)
    return target_quaternion(args.static_args.gs_lat, args.static_args.gs_lon, args.dynamic_args)
end
function mode_quaternion(::Type{CamPointing}, args::PointingArguments)
    return target_quaternion(args.static_args.cam_lat, args.static_args.cam_lon, args.dynamic_args)
end

function target_quaternion(lat, lon, args)
    lat = deg2rad(lat)
    lon = deg2rad(lon)
    # WGS-84 ellipsoidal model parameters
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
    return align_frame_with_vector(vec_cubesat_to_location_body, normalize(args.sun_body), [1.0,0,0], [0,0,-1.0] )
end
