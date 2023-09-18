abstract type AbstractSensor end

@kwdef @concrete struct NadirSensor <: AbstractSensor
    vfov = 0.6981317007977318
    hfov = 0.7853981633974483
    maximum_rate = 0.24434609527920614
    σ_cross_sight = 0.005817764173314432
    σ_roll = 0.005817764173314432
    q_body_sensor::Quaternion
end

@kwdef @concrete struct StarTracker <: AbstractSensor
    fov = 0.7330382858376184
    maximum_rate = 0.005235987755982988
    σ_cross_sight = 0.00011635528346628864
    σ_roll = 0.00034906585039886593
    q_body_sensor::Quaternion
end

@kwdef @concrete struct SunSensor <: AbstractSensor
    fov = 1.4835298641951802
    maximum_rate = 1.2217304763960306
    σ_cross_sight = 0.0017453292519943296
    σ_roll = 0.0017453292519943296
    q_body_sensor::Quaternion
end
