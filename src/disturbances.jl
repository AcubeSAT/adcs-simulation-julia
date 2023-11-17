"""
======================================================================== 
   This function calculates the residual magnetic torque
 
   Inputs:
     cosines      - Cosine of angle between body frame axes and orbit frame z-axis
     B_body       - Magnetic field expressed in body frame

   Ouputs:
     res_mag_tor  - Residual magnetic torque
     rm           - Residual magnetic moment
 
 ======================================================================== 
"""

function residual_magnetic_torque(cosines, B_body)
    mean_rm_base = [0.05, 0.05, 0.05]
    sign_vector = sign.(-cosines)
    rm = mean_rm_base + 0.1 * mean_rm_base .* sign_vector + 0.0005 * rand(3)
    res_mag_tor = cross(rm, B_body)

    return res_mag_tor, rm
end


"""
 ======================================================================== 
   This function calculates the solar pressure that acts upon
   the spacecraft.
 
   Inputs:
     R_OB          - Transformation matrix from orbit to body frame
     sun_vector    - Unit vector from satellite to sun expressed in orbit frame
     cm            - Center of mass
 
   Ouputs:
     solar_rad_tor - Solar pressure torque
     area          - Satellite's area projected to the sun
     cosines       - Cosine of the angles between the body frame axes of the satellite and the sun vector transformed into the body frame

 ======================================================================== 
"""

function solar_radiation_pressure(R_OB, sun_vector, cm)
    Fs = 1367;                                          "Solar Constant [W/m^2]"
    c  = 3e8;                                           "Speed of light [m/s]"
    reflectance_factor = 0.6;

    y_0 = -R_OB * sun_vector;                           "Sun vector in body frame"
    y_0 = y_0 / norm(y_0);                              "Make it unit vector"

    proj_Xb_Zo = dot([1 0 0], y_0) * y_0;               "Projection of each body frame axes unit vector to orbit frame z axis" 
    proj_Yb_Zo = dot([0 1 0], y_0) * y_0;           
    proj_Zb_Zo = dot([0 0 1], y_0) * y_0;

    proj_Xb_XYo = [1, 0, 0] - proj_Xb_Zo;               "Projection of each body frame axes unit vector to orbit frame x,y plane"
    proj_Yb_XYo = [0, 1, 0] - proj_Yb_Zo;
    proj_Zb_XYo = [0, 0, 1] - proj_Zb_Zo;

    Ax = 0.34 * norm(cross(proj_Yb_XYo, proj_Zb_XYo));  "Surface projections to orbit frame x,y plane"
    Ay = 0.34 * norm(cross(proj_Xb_XYo, proj_Zb_XYo));
    Az = 0.1 * norm(cross(proj_Xb_XYo, proj_Yb_XYo));  
    area = [Ax Ay Az];
    
    cos_Xb_Xo = [1 0 0] * y_0;
    cos_Yb_Yo = [0 1 0] * y_0;
    cos_Zb_Zo = [0 0 1] * y_0;
    cosines = [cos_Xb_Xo, cos_Yb_Yo, cos_Zb_Zo]; 

    solar_pressure_center = Diagonal([0.05, 0.05, 0.17]); 

    solar_pressure_center[1,:] = sign.(-cos_Xb_Xo) .* solar_pressure_center[1,:]; 
    solar_pressure_center[2,:] = sign.(-cos_Yb_Yo) .* solar_pressure_center[2,:]; 
    solar_pressure_center[3,:] = sign.(-cos_Zb_Zo) .* solar_pressure_center[3,:]; 

    T1 = (Fs / c) * Ax * (1 + reflectance_factor) * cross(y_0, solar_pressure_center[1,:] - cm);
    T2 = (Fs / c) * Ay * (1 + reflectance_factor) * cross(y_0, solar_pressure_center[2,:] - cm);
    T3 = (Fs / c) * Az * (1 + reflectance_factor) * cross(y_0, solar_pressure_center[3,:] - cm);

    solar_rad_tor = (T1+T2+T3);

    return solar_rad_tor, area, cosines
end