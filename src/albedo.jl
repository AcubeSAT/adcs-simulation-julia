@kwdef @concrete mutable struct AlbedoParameters
    satECEF                    # ECEF position vector of satellite
    sunECEF                    # ECEF position vector of Sun
    TOMSMatrix                 # Matrix containing the reflectivity value of each grid, from NASA TOMS project
    radius = 6371.01 * 10^3    # Mean Earth radius in meters
    irr = 1                    # Solar irradiance constant  
    TOMSrows = 180             # TOMS' matrix number of rows
    TOMScolumns = 288          # TOMS' matrix number of columns
    dx = 2π / TOMScolumns
    dy = π / TOMSrows
    dPhi = deg2rad(180 / TOMSrows)
    dTheta = deg2rad(360 / TOMScolumns)
    dPhiHalf = dPhi / 2 
end

function rad2ind(theta, phi, AP::AlbedoParameters)
    i = round((π - phi - AP.dy / 2) / AP.dy)
    j = round((π + theta - AP.dx / 2) / AP.dx)

    k = firstindex(AP.TOMSMatrix)
    return (ifelse(i < k, k, i), ifelse(j < k, k, j))
end

function ind2rad(i, j, AP::AlbedoParameters)
    return (theta, phi) = (-π + AP.dx / 2 + j * AP.dx, π - AP.dy / 2 - i * AP.dy)
end

function cellArea(i, j, AP::AlbedoParameters)
    radians = ind2rad(i, j, AP)

    maxPhi = radians[2] + AP.dPhiHalf
    minPhi = radians[2] - AP.dPhiHalf

    return AP.radius * AP.radius * AP.dTheta * (cos(minPhi) - cos(maxPhi))
end 

function gridAngle(i, j, iSun, jSun, AP::AlbedoParameters)
    loopRad = ind2rad(i, j, AP)
    sunRad = ind2rad(iSun, jSun, AP)

    return acos(sin(loopRad[2]) * sin(sunRad[2]) * cos(loopRad[1] - sunRad[1]) + cos(loopRad[2]) * cos(sunRad[2]))
end 

function calculateAlbedo(AP::AlbedoParameters)
    sunSph = SphericalFromCartesian()(AP.sunECEF)
    sunECEFSpherical = [sunSph.r, sunSph.θ, sunSph.ϕ]
    sunECEFSpherical[2] = π / 2 - sunECEFSpherical[2]
    indSun = rad2ind(sunECEFSpherical[1], sunECEFSpherical[2], AP)

    albedo = zeros(AP.TOMSrows, AP.TOMScolumns)

    for i in 1 : AP.TOMSrows
        for j in 1 : AP.TOMScolumns

            # calculate angle of incident irradiance
            irrAngle = gridAngle(i, j, indSun[1], indSun[2], AP)       
            irrAngle = min(irrAngle, π / 2)
            
            # calculate incident power 
            power = AP.irr * cellArea(i, j, AP) * cos(irrAngle)

            # calculate cartesian position of grid in ECEF
            gridRad = ind2rad(i, j, AP)
            gridTheta = gridRad[1]
            gridPhi = gridRad[2]
            grid = CartesianFromSpherical()(Spherical(AP.radius, gridTheta, π / 2 - gridPhi))

            # calculate the distance from grid to satellite  
            satDist = norm(AP.satECEF - grid)
            
            # calculate angle between satellite and grid
            satGridAngle = acos(dot(AP.satECEF - grid, grid) / (satDist * norm(grid)))

            # finally, calculate albedo value 
            albedo[i, j] = power * AP.TOMSMatrix[i, j] * cos(satGridAngle) / (π * satDist^2)

        end
    end

    return albedo
end

# Earth Albedo Model based on
#       Bhanderi, D. D. V. (2005). 
#       Spacecraft Attitude Determination with Earth Albedo Corrected Sun Sensor Measurements. 
#       Aalborg University.