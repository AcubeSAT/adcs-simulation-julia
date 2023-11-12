module ADCSSims

using Accessors
using ConcreteStructs

using StaticArrays

using SatelliteDynamics
using SatelliteToolboxGeomagneticField

using Infiltrator
using LinearAlgebra
using Plots
using Random
using Statistics

include("quaternion.jl")
include("sensor.jl") # FIXME: This is not needed for AcubeSAT, given its placement
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

include("disturbances.jl")

export Quaternion, QuaternionF16, QuaternionF32, QuaternionF64, rotvec, from_rotation_matrix
export NadirSensor, StarTracker, SunSensor, in_fov, available, qerr, emulate_estimation
export ReactionWheel, stribeck, deadzone_compensation, saturation_compensation
export PDController, calculate_torque, decompose_torque
export mse, qloss

end
