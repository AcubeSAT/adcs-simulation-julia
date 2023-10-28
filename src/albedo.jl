@kwdef @concrete struct AlbedoParameters
    satPosition # ECEF position vector of satellite
    sunPosition # ECEF position vector of Sun
    reflectivityData # reflectivity value of each grid
    EMR = 6371.01 * 10^3 # Earth Mean Radius in meters
    TOMSrows = 180
    TOMScolumns = 288
    dx = 2π / TOMScolumns
    dy = π / TOMSrows
    dPhi = deg2rad(180 / TOMSrows)
    dTheta = deg2rad(360 / TOMScolumns)
    dPhiHalf = dPhi / 2 
    solarIrradiance = 1
end

function rad2ind(theta, phi, AP::AlbedoParameters)
    i = round((π - phi - AP.dy / 2) / AP.dy)
    j = round((π + theta - AP.dx / 2) / AP.dx)

    (i, j) = (max(0, i), max(0, j))
end

function ind2rad(i, j, AP::AlbedoParameters)
    (theta, phi) = (-π + AP.dx / 2 + j * AP.dx, π - AP.dy / 2 - i * AP.dy)
end

function calculateCellArea(i, j, AP::AlbedoParameters)
    radians = ind2rad(i, j, AP)

    maxPhi = radians[2] + AP.dPhiHalf
    minPhi = radians[2] - AP.dPhiHalf

    area = AP.EMR * AP.EMR * deltaTheta * (cos(minPhi) - cos(maxPhi))
end 

function gridAngle(loopI, loopJ, sunIndexI, sunIndexJ, AP::AlbedoParameters)
    loopRadians = ind2rad(loopI, loopJ, AP)
    sunRadians = ind2rad(sunIndexI, sunIndexJ, AP)

    angle = acos(sin(loopRadians[2]) * sin(sunRadians[2])* cos(loopRadians[1] - sunRadians[1]) + cos(loopRadians[2]) * cos(sunRadians[2]))
end 

# All vector parameters of the following function must be expressed in ECEF frame
function calculateAlbedo(AP::AlbedoParameters)
    sunPositionSpherical = SphericalFromCartesian()(AP.sunPosition)
    sunPositionSpherical[2] = π / 2 - sunPositionSpherical[2]

    sunIndices = rad2ind(sunPositionSpherical[1], sunPositionSpherical[2], AP)

    for p in 1:AP.TOMSrows
        for k in 1:AP.TOMScolumns

            angleOfIncidentSolarIrradiance = gridAngle(p, k, sunIndices[1], sunIndices[2], AP)

            angleOfIncidentSolarIrradiance = min(angleOfIncidentSolarIrradiance, π / 2 )
            
            incidentPower = AP.solarIrradiance * calculateCellArea(p, k, AP) * cos(angleOfIncidentSolarIrradiance)
            gridRadians = ind2rad(p, k, AP)
            gridTheta = gridRadians[1]
            gridPhi = gridRadians[2]

            grid = CartesianFromSpherical()([AP.EMR, gridTheta, π / 2 - gridPhi])

            satelliteDistance = norm(AP.satPosition - grid)
            
            satelliteGridAngle = acos((transpose((AP.satPosition - grid) / satelliteDistance)) * grid / norm(grid))

            albedo[p, k] = incidentPower * AP.refectivityData[p, k] * cos(satelliteGridAngle) / (π * satelliteDistance^2)

        end
    end

    return albedo
end

# Earth Albedo Model based on
#       Bhanderi, D. D. V. (2005). 
#       Spacecraft Attitude Determination with Earth Albedo Corrected Sun Sensor Measurements. 
#       Aalborg University.