module ADCSSims

using ConcreteStructs

using StaticArrays

using SatelliteDynamics
using SatelliteToolboxGeomagneticField

using LinearAlgebra
using Plots
using Statistics
using CoordinateTransformations
using PoissonRandom

include("quaternion.jl")
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
include("albedo.jl")

export Quaternion, QuaternionF16, QuaternionF32, QuaternionF64, rotvec
export ReactionWheel, stribeck, deadzone_compensation, saturation_compensation
export PDController, calculate_torque
export mse, qloss

end
