abstract type AbstractSensor end

@kwdef @concrete struct NadirSensor <: AbstractSensor
    vfov = 0.6981317007977318
    hfov = 0.7853981633974483
    maximum_rate = 0.24434609527920614
    σ_cross_sight = 0.005817764173314432
    σ_roll = 0.005817764173314432
    q_body_sensor = Quaternion(sqrt(2)/2, 0.0, sqrt(2)/2, 0.0)
    abs_weight = 0.333333
end

@kwdef @concrete struct StarTracker <: AbstractSensor
    fov = 0.7330382858376184
    maximum_rate = 0.005235987755982988
    σ_cross_sight = 0.00011635528346628864
    σ_roll = 0.00034906585039886593
    q_body_sensor = Quaternion(sqrt(2)/2, 0.0, sqrt(2)/2, 0.0)
    abs_weight = 227.272727
end

@kwdef @concrete struct SunSensor <: AbstractSensor
    fov = 1.4835298641951802
    maximum_rate = 1.2217304763960306
    σ_cross_sight = 0.0017453292519943296
    σ_roll = 0.0017453292519943296
    q_body_sensor = QuaternionF64(1,0,0,0)
    abs_weight = 8.333333
end

# FOV should be in radians, half of the sensor FOV
function in_fov(tpos, fov)
    LinearAlgebra.normalize!(tpos)
    return acos(dot(tpos, [0 0 1])) <= fov
end

function in_fov(tpos, vfov, hfov)
    LinearAlgebra.normalize!(tpos)
    θv = abs(acos(tpos[3]))
    θh = abs(atan(tpos[2], tpos[1]))
    return θv <= vfov && θh <= hfov
end

# function available(NS::NadirSensor, target_vector, w)
#     return w <= NS.maximum_rate && in_fov(target_vector, NS.vfov, NS.hfov)
# end

# function available(ST::StarTracker, target_vector, w)
#     return w <= ST.maximum_rate && !in_fov(target_vector, ST.fov)
# end

# function available(SN::SunSensor, target_vector, w)
#     return w <= SN.maximum_rate && in_fov(target_vector, SN.fov)
# end

function available(NS::NadirSensor, target_vector, w)
    return in_fov(target_vector, NS.vfov, NS.hfov)
end

function available(ST::StarTracker, target_vector, w)
    return !in_fov(target_vector, ST.fov)
end

function available(SN::SunSensor, target_vector, w)
    return in_fov(target_vector, SN.fov)
end

available(::AbstractSensor) = error("available is not defined in the abstract type")
function qerr(S::AbstractSensor)
    δθx = randn() * S.σ_cross_sight
    δθy = randn() * S.σ_cross_sight
    δθz = randn() * S.σ_roll
    return LinearAlgebra.normalize(Quaternion(1.0, δθx / 2, δθy / 2, δθz / 2))
end

function emulate_estimation(sensors, target_vectors, w)
    err_qs = Quaternion[]
    abs_weights = typeof(sensors[1].abs_weight)[]
    for (sensor, tvec) in zip(sensors, target_vectors)
        if available(sensor, tvec, norm(w))
            push!(err_qs, qerr(sensor))
            push!(abs_weights, sensor.abs_weight)
        end
    end
    @assert !isempty(err_qs)
    rel_weights = abs_weights ./ sum(abs_weights)
    return true, normalize(sum(err_qs .* rel_weights))
end
