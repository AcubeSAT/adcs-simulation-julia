module ADCSSims

using Accessors
using ConcreteStructs

using StaticArrays

using CoordinateTransformations

using SatelliteDynamics
using SatelliteToolboxGeomagneticField
using SatelliteToolboxGravityModels

using BenchmarkTools
using Infiltrator

using Dates
using LinearAlgebra
using Plots
using Random
using Statistics

include("quaternion.jl")
include("pointing_modes.jl")
include("simulation.jl")
include("control.jl")
include("disturbances.jl")
include("sensor.jl")
include("reaction_wheel.jl")
include("dynamics.jl")
include("losses.jl")
include("utils.jl")
include("parameters.jl")
include("plots.jl")
include("frames.jl")

export Quaternion,
    QuaternionF16,
    QuaternionF32,
    QuaternionF64,
    rotvec,
    from_rotation_matrix,
    to_rotation_matrix,
    to_euler_angles,
    align_frame_with_vector
export NadirSensor,
    StarTracker,
    SunSensor,
    in_fov,
    available,
    qerr,
    emulate_estimation
export PointingMode, NadirPointing, GroundTargetPointing, SunPointing, PointingArguments, mode_quaternion
export ReactionWheel, stribeck, deadzone_compensation, saturation_compensation
export PDController, calculate_torque, decompose_torque
export mse, qloss
export SimulationParams, SimulationContext
export parse_pointing_mode
end
