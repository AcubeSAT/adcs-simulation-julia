include("parameters.jl")
include("orbitFunctions.jl")

timestep = 0.1
totalTime = 8*5545.0
timeOffsets = 0:timestep:totalTime

propagatorsArray = orbitInit("TLE.txt",size(timeOffsets)[1]) 

magneticFieldECI, magneticFieldOrbit, sunPositionECI, sunPositionOrbit, eclipseFlag =
    getReferenceVectors(propagatorsArray, timeOffsets);

#TODO parameters + documentation