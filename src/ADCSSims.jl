module ADCSSims

using SatelliteDynamics
using SatelliteToolboxGeomagneticField

using LinearAlgebra
using Plots
using Statistics

include("dynamics.jl")
include("utils.jl")
include("parameters.jl")
include("mekf.jl")
include("noise_models.jl")
include("simulation.jl")
include("plots.jl")

end