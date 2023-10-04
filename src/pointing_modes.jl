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