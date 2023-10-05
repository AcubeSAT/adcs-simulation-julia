abstract type PointingMode end
struct SunPointing <: PointingMode end
struct NadirPointing <: PointingMode end
struct GSPointing <: PointingMode end

@concrete struct PointingArguments
    sun_body
    nadir_body
    qeci2body
    qeci2orbit
end

function mode_quaternion(::Type{PointingMode}, args::PointingArguments)
    error("Mode not implemented")
end
function mode_quaternion(::Type{SunPointing}, args::PointingArguments)
    return align_frame_with_vector(args.sun_body, args.nadir_body, [0,0,-1], [0,1,0])
end
function mode_quaternion(::Type{NadirPointing}, args::PointingArguments)
    return args.qeci2body * conj(args.qeci2orbit)
end