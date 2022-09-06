include("utilities.jl")

function NED2ECEF(vectorNED, latitude::Float64, longtitude::Float64)
    c1 = [-sin(latitude) * cos(longtitude); -sin(latitude) * sin(longtitude); cos(latitude)]
    c2 = [-sin(longtitude); cos(longtitude); 0]
    c3 = [-cos(latitude) * cos(longtitude); -cos(latitude) * sin(longtitude); -sin(latitude)]
    R = hcat(c1, c2, c3)
    q = DCM2Quat(R)
    vectorECEF = rotateVector(q, vectorNED)
end

function ECI2ECEFWrapper(julianDate::Number)
    r_eci_to_ecef(DCM, GCRF(), ITRF(), julianDate, eop_IAU2000A)
end

function ECEF2ECI(vectorECEF::Vector{Float64}, timeGST::Float64)
    CGAST = cos(-timeGST)
    SGAST = sin(-timeGST)
    R = vcat([CGAST SGAST 0], [-SGAST CGAST 0], [0 0 1])
    q = DCM2Quat(R)
    vectorECI = rotateVector(q, vectorECEF)
end

function ECI2Orbit(vectorECI::Vector{Float64}, propagator::OrbitPropagatorSGP4{Float64})

    inclination = propagator.sgp4d.i_k
    rightAscensionOfAscendingNode = propagator.sgp4d.Ω_k
    argumentOfPerigee = propagator.sgp4d.ω_k

    R = zeros(3, 3)
    R[1, 1] =
        -sin(argumentOfPerigee) * sin(rightAscensionOfAscendingNode) * cos(inclination) +
        cos(argumentOfPerigee) * cos(rightAscensionOfAscendingNode)
    R[1, 2] =
        sin(argumentOfPerigee) * cos(inclination) * cos(rightAscensionOfAscendingNode) +
        sin(rightAscensionOfAscendingNode) * cos(argumentOfPerigee)
    R[1, 3] = sin(inclination) * sin(argumentOfPerigee)

    R[2, 1] =
        -sin(argumentOfPerigee) * cos(rightAscensionOfAscendingNode) -
        sin(rightAscensionOfAscendingNode) * cos(inclination) * cos(argumentOfPerigee)
    R[2, 2] =
        -sin(argumentOfPerigee) * sin(rightAscensionOfAscendingNode) +
        cos(inclination) * cos(argumentOfPerigee) * cos(rightAscensionOfAscendingNode)
    R[2, 3] = sin(inclination) * cos(argumentOfPerigee)

    R[3, 1] = sin(inclination) * sin(rightAscensionOfAscendingNode)
    R[3, 2] = -sin(inclination) * cos(rightAscensionOfAscendingNode)
    R[3, 3] = cos(inclination)

    vectorOrbit = R * vectorECI
end
