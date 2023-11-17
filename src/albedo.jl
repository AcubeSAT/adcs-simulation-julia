@kwdef @concrete mutable struct AlbedoParameters
    sat_ecef                    # ECEF position vector of satellite
    sun_ecef                    # ECEF position vector of Sun
    toms_matrix                 # Matrix containing the reflectivity value of each grid, from NASA TOMS project
    radius = 6371.01 * 10^3     # Mean Earth radius in meters
    irr = 1                     # Solar irradiance constant  
    toms_rows = 180             # TOMS' matrix number of rows
    toms_columns = 288          # TOMS' matrix number of columns
    dx = 2π / toms_columns
    dy = π / toms_rows
    dphi = deg2rad(180 / toms_rows)
    dtheta = deg2rad(360 / toms_columns)
    dphi_half = dphi / 2 
end

function rad2ind(theta, phi, AP::AlbedoParameters)
    i = round((π - phi - AP.dy / 2) / AP.dy)
    j = round((π + theta - AP.dx / 2) / AP.dx)

    k = firstindex(AP.toms_matrix)
    return ifelse(i < k, k, i), ifelse(j < k, k, j)
end

function ind2rad(i, j, AP::AlbedoParameters)
    return theta, phi = -π + AP.dx / 2 + j * AP.dx, π - AP.dy / 2 - i * AP.dy
end

function cell_area(i, j, AP::AlbedoParameters)
    radians = ind2rad(i, j, AP)

    max_phi = radians[2] + AP.dphi_half
    min_phi = radians[2] - AP.dphi_half

    return AP.radius^2 * AP.dtheta * (cos(min_phi) - cos(max_phi))
end 

function grid_angle(i, j, i_sun, j_sun, AP::AlbedoParameters)
    loop_rad = ind2rad(i, j, AP)
    sun_rad = ind2rad(i_sun, j_sun, AP)

    return acos(sin(loop_rad[2]) * sin(sun_rad[2]) * cos(loop_rad[1] - sun_rad[1]) + cos(loop_rad[2]) * cos(sun_rad[2]))
end 

function calculate_albedo(AP::AlbedoParameters)
    sun_sph = SphericalFromCartesian()(AP.sun_ecef)
    sun_ecef_sph = [sun_sph.r, sun_sph.θ, sun_sph.ϕ]
    sun_ecef_sph[2] = π / 2 - sun_ecef_sph[2]
    ind_sun = rad2ind(sun_ecef_sph[1], sun_ecef_sph[2], AP)

    albedo = zeros(AP.toms_rows, AP.toms_columns)

    for i in 1 : AP.toms_rows
        for j in 1 : AP.toms_columns

            # calculate angle of incident irradiance
            irr_angle = grid_angle(i, j, ind_sun[1], ind_sun[2], AP)       
            irr_angle = min(irr_angle, π / 2)
            
            # calculate incident power 
            power = AP.irr * cell_area(i, j, AP) * cos(irr_angle)

            # calculate cartesian position of grid in ECEF
            grid_rad = ind2rad(i, j, AP)
            grid_theta = grid_rad[1]
            grid_phi = grid_rad[2]
            grid = CartesianFromSpherical()(Spherical(AP.radius, grid_theta, π / 2 - grid_phi))

            # calculate the distance from grid to satellite  
            sat_dist = norm(AP.sat_ecef - grid)
            
            # calculate angle between satellite and grid
            sat_grid_angle = acos(dot(AP.sat_ecef - grid, grid) / (sat_dist * norm(grid)))

            # finally, calculate albedo value 
            albedo[i, j] = power * AP.toms_matrix[i, j] * cos(sat_grid_angle) / (π * sat_dist^2)

        end
    end

    return albedo
end

# Earth Albedo Model based on
#       Bhanderi, D. D. V. (2005). 
#       Spacecraft Attitude Determination with Earth Albedo Corrected Sun Sensor Measurements. 
#       Aalborg University.