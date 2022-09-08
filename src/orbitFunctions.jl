using SatelliteToolbox
using Quaternions
using Plots
using JLD2
using Infiltrator

include("sunFunctions.jl")
include("transformations.jl")
include("utilities.jl")

const eop_IAU2000A = load_object("IGRF_data/eop_IAU2000A.jld2")
const monthsInYear = 12.0
const averageMonthDays = 30.436875
const hoursInDay = 24.0
const minutesInHour = 60.0
const secondsInMinute = 60.0
const secondsInDay = 86400.0

"""
Initialize an orbit propagator object given the path to a TLE file
# Arguments
- `tleString::String`: path to TLE file.
"""
function orbitInit(tleString::String, copies=1)
    tle = read_tle(tleString)[1]
    propagatorsArray = [init_orbit_propagator(Val(:sgp4), tle) for copy in range(1, copies)]
end

"""
Propagate the orbit and return the magnetic field, the sun vectors and eclipse flags
# Arguments
- 'propagator': propagator object created by orbitInit
- 'timestep': propagation timestep
- 'totalTime': length of time for which to run
"""
function getReferenceVectors(
    propagatorsArray::Vector{OrbitPropagatorSGP4{Float64}},
    timeOffsets,
)

    initialJulianDate = propagatorsArray[1].sgp4d.epoch

    julianDate = initialJulianDate .+ (timeOffsets / secondsInDay)

    gregorianDate = jd_to_date.(julianDate)

    years = [x[1] for x in gregorianDate]
    months = [x[2] for x in gregorianDate]
    days = [x[3] for x in gregorianDate]

    decimalGregorianDate =
        years .+ (months[2] .- 1) / monthsInYear .+
        (days .- 1) / averageMonthDays / monthsInYear


    satelliteDataECI = propagate!.(propagatorsArray, timeOffsets)
    positionECI = [Vector{Float64}(x[1]) for x in satelliteDataECI]

    rotationECI2ECEF = ECI2ECEFWrapper.(julianDate)

    positionECEF = rotationECI2ECEF .* positionECI
    positionLLA = ecef_to_geodetic.(positionECEF)

    latitude = [x[1] for x in positionLLA]
    longtitude = [x[2] for x in positionLLA]
    altitude = [x[3] for x in positionLLA]

    magneticFieldNED =
        igrf.(decimalGregorianDate, altitude, latitude, longtitude, Val(:geodetic))
    magneticFieldECEF = ned_to_ecef.(magneticFieldNED, latitude, longtitude, altitude)
    magneticFieldECEF = Vector{Vector}(magneticFieldECEF)

    greenwichSiderialTime = julianDay2GST.(julianDate)
    magneticFieldECI = ECEF2ECI.(magneticFieldECEF, greenwichSiderialTime)
    magneticFieldOrbit = ECI2Orbit.(magneticFieldECI, propagatorsArray)

    sunPositionECI = sunPosition.(julianDate)

    sunPositionOrbit = ECI2Orbit.(sunPositionECI, propagatorsArray)

    eclipseFlag = calculateEclipse.(positionECI, sunPositionECI)

    returnVectors = (
        magneticFieldECI,
        magneticFieldOrbit,
        sunPositionECI,
        sunPositionOrbit,
        eclipseFlag,
    )

end
