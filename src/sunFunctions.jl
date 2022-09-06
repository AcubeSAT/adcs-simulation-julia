using LinearAlgebra

"""
Calculate sun position given time
# Arguments
- `julianDay::Float64`: time in JD format.
"""
function sunPosition(julianDay::Float64)

    time = julianDay
    # Convert time
    ut1 = (time - 2451545.0) / 36525.0

    # ------ Algorithm --------
    meanlong = 280.4606184 + 36000.77005361 * ut1
    meanlong = mod(meanlong, 360) #deg
    meananomaly = 357.5277233 + 35999.05034 * ut1
    meananomaly = mod(meananomaly * pi / 180, 2 * pi) #rad

    if meananomaly < 0
        meananomaly = 2 * pi + meananomaly
    end

    eclplong = meanlong + 1.91466471 * sin(meananomaly) + 0.019994643 * sin(2 * meananomaly) #deg
    obliquity = 23.439291 - 0.0130042 * ut1 #deg
    meanlong = meanlong * pi / 180 #rad

    if meanlong < 0
        meanlong = 2 * pi + meanlong
    end

    eclplong = eclplong * pi / 180 #rad
    obliquity = obliquity * pi / 180 #rad

    # Magnitude of sun vector and ECI calculations
    magr = 1.000140612 - 0.016708617 * cos(meananomaly) - 0.000139589 * cos(2 * meananomaly) #au's
    sunPositionECI = [magr * cos(eclplong); magr * cos(obliquity) * sin(eclplong); magr * sin(obliquity) * sin(eclplong)]

end

"""
Calculate whether the satellite is in eclipse
# Arguments
- `satellitePositionECI::Float64`
- `sunPositionECI::Float64`

# Returns
- `Eclipse flag:: Boolean`

"""
function calculateEclipse(satellitePositionECI::Vector{Float64}, sunPositionECI::Vector{Float64})

    earthRadius = 6371
    sunRadius = 696000
    AU = 149600000

    x1 = earthRadius * AU / (sunRadius + earthRadius)
    alpha1 = pi - acos(earthRadius / x1) - acos(earthRadius / norm(satellitePositionECI))

    x2 = earthRadius * AU / (sunRadius - earthRadius)
    alpha2 = acos(earthRadius / x2) - acos(earthRadius / norm(satellitePositionECI))

    alpha = pi - acos(dot(sunPositionECI, satellitePositionECI) / (norm(sunPositionECI) * norm(satellitePositionECI)))

    if alpha2 < alpha && alpha < alpha1
        eclipse = true
    elseif alpha < alpha2
        eclipse = true
    else
        eclipse = false
    end

    eclipse
end