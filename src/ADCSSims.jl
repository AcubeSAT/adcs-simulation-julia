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

include("pointing_modes.jl")
include("quaternion.jl")
include("disturbances.jl")
include("sensor.jl")
include("reaction_wheel.jl")
include("control.jl")
include("noise_models.jl")
include("dynamics.jl")
include("losses.jl")
include("utils.jl")
include("parameters.jl")
include("mekf.jl")
include("simulation.jl")
include("plots.jl")
include("frames.jl")

export Quaternion, QuaternionF16, QuaternionF32, QuaternionF64, rotvec, from_rotation_matrix, to_rotation_matrix, to_euler_angles, align_frame_with_vector
export NadirSensor, StarTracker, SunSensor, in_fov, available, qerr, emulate_estimation
export PointingMode, NadirPointing, GSPointing, SunPointing,CamPointing, PointingArguments, mode_quaternion
export ReactionWheel, stribeck, deadzone_compensation, saturation_compensation
export PDController, calculate_torque, decompose_torque
export mse, qloss
export SimulationParams, SimulationContext
end
