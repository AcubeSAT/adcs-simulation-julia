module ADCSSims

using Accessors
using ConcreteStructs

using StaticArrays

using SatelliteDynamics
using SatelliteToolboxGeomagneticField

using LinearAlgebra
using Plots
using Random
using Statistics

include("quaternion.jl")
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

export Quaternion
export NadirSensor, StarTracker, SunSensor, in_fov, available, qerr, estimateq
export ReactionWheel, stribeck, deadzone_compensation, saturation_compensation
export PDController, calculate_torque, decompose_torque
export mse, qloss

end
