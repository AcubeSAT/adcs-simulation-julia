include("parameters.jl")
include("orbitFunctions.jl")

propagator = orbitInit("TLE.txt")
magneticFieldECI, magneticFieldOrbit, sunPositionECI, sunPositionOrbit, eclipseFlag =
    getReferenceVectors(propagator, 0.1, 5545.0);

#TODO parameters + documentation