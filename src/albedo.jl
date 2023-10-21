earthRadius = 6371.01 * 10^3
TOMSrows = 180
TOMScolumns = 288
dx = 2π / TOMScolumns
dy = π / TOMSrows
solarIrradiance = 1

function rad2ind(theta, phi)
    i = round((π - phi - dy / 2) / dy)
    j = round((π + theta - dx / 2) / dx)

    (i, j) = (max(0, i), max(0, j))
end

function ind2rad(index_i, index_j)
    (theta, phi) = (-π + dx / 2 + index_j * dx, π - dy / 2 - index_i * dy)
end

function calculateCellArea(index_i, index_j)
    radians = ind2rad(index_i, index_j)

    deltaPhi = deg2rad(180 / TOMSrows)
    deltaTheta = deg2rad(360 / TOMScolumns)

    maxPhi = radians[2] + deltaPhi / 2
    minPhi = radians[2] - deltaPhi / 2

    area = earthRadius * earthRadius * deltaTheta * (cos(minPhi) - cos(maxPhi))
end 

function gridAngle(loopI, loopJ, sunIndex_i, sunIndex_j)
    loopRadians = ind2rad(loopI, loopJ)
    sunRadians = ind2rad(sunIndex_i, sunIndex_j)

    angle = acos(sin(loopRadians[2]) * sin(sunRadians[2])* cos(loopRadians[1] - sunRadians[1]) + cos(loopRadians[2]) * cos(sunRadians[2]))
end 

# All vector parameters of the following function must be expressed in ECEF frame
function calculateAlbedo(satPosition, sunPosition, refectivityData)
    sunPositionSpherical = SphericalFromCartesian()(sunPosition)
    sunPositionSpherical[2] = π / 2 - sunPositionSpherical[2]

    sunIndices = rad2ind(sunPositionSpherical[1], sunPositionSpherical[2])

    for i in 1:TOMSrows
        for j in 1:TOMScolumns

            angleOfIncidentSolarIrradiance = gridAngle(i, j, sunIndices[1], sunIndices[2])

            angleOfIncidentSolarIrradiance = min(angleOfIncidentSolarIrradiance, π / 2 )
            
            incidentPower = solarIrradiance * calculateCellArea(i, j) * cos(angleOfIncidentSolarIrradiance)
            gridRadians = ind2rad(i, j)
            gridTheta = gridRadians[1]
            gridPhi = gridRadians[2]

            grid = CartesianFromSpherical()([gridTheta, π / 2 - gridPhi, earthRadius])

            satelliteDistance = norm(satPosition - grid)
            
            satelliteGridAngle = acos((transpose((satPosition - grid) / satelliteDistance)) * grid / norm(grid))

            albedo[i, j] = incidentPower * refectivityData[i, j] * cos(satelliteGridAngle) / (π * satelliteDistance^2)

        end
    end

    return albedo
end

# Earth Albedo Model based on
#       Bhanderi, D. D. V. (2005). 
#       Spacecraft Attitude Determination with Earth Albedo Corrected Sun Sensor Measurements. 
#       Aalborg University.