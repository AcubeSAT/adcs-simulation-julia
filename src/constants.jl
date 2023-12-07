#=====================================================================
In this function all parameters regarding the physical architecture 
of the satellite and the ADCS components are set.
    Input: Nothing
    Output: Struct where all values of the parameters are included

    Parameters:
        m             - Satellite mass [kg]
        lx            - Satellite X-axis length [m]
        ly            - Satellite Y-axis length [m]
        lz            - Satellite Z-axis length [m]
        N_coils       - Number of turns of the coils
        A_coils       - Total area of the coils [m^2]
        R_coils       - Coils resistance [Ohm]
        mt            - Proportional parameter for calculating the current on each MTQ
                        For torquer rods: mt = 1.66 * (length/diameter)^1.5.
                        For air core magnetorquers: mt = 1.
        W_up_limit    - Upper limit of absolute power of mtq's [W]
        I_up_limit    - Upper limit of absolute current of mtq's [A]
        Km            - Motor torque constant (RW)
        Rb            - Armature resistance (RW)
        Kv            - Motor velocity constant (back-EMF constant) -> 1/Km
        b_friction    - Viscuous friction
        c_friction    - Coulomb friction
        Jw            - Inertia of the RW
        A             - Parameter which indicates the percentance of RW desaturation torque
        lim_dz        - The absolute value of the limits of the deadzone (in rpm)
        p_400         - Atmospheric density in an altitude of 400km [kg/m^3]
        p_500         - Atmospheric density in an altitude of 500km [kg/m^3]
        Ix            - Xaxis inertia (simplified)
        Iy            - Yaxis inertia (simplified)
        Iz            - Zaxis inertia (simplified)
        PMI           - Principal Moments of Inertia
        PAI           - Principal Axes of Inertia
        Cm            - Center of mass
        I             - Inertia matrix -> PAI*PMI*PAI'
        I_inv 	      - Inverse inertia matrix -> eye(3,3)/I;
        Re            - Earth radius [m]
        Rs            - Satellite altitude [m]
        Radius        - Distance from earth center to satellite [m]
        G             - Earth gravitational constant
        M             - Earth mass
        w_o           - Satellite angular velocity relative to Earth
        v_satellite   - Satellite's velocity in orbit
        known_rm      - Residual magnetic dipole
        orbitPeriod   - Orbit Period
        n             - Total simulation time
        mtq_max       - Maximum magnetic dipole of magnetorquers
        rw_max_torque - Maximum torque of the Reaction Wheel

=====================================================================#

struct SatelliteParameters
    Radius::Float64
    w_o::Float64
    v_satellite::Float64
    n::Float64
    orbitPeriod::Float64
    I::Matrix{Float64}
    I_inv::Matrix{Float64}
    w_o_io::Vector{Float64}        # what is this for????
    N_coils::Vector{Int}
    A_coils::Vector{Float64}
    R_coils::Vector{Float64}
    mt::Vector{Float64}
    I_up_limit::Vector{Float64}
    W_up_limit::Vector{Float64}
    Jw::Float64
    A::Float64
    Km::Float64
    Kv::Float64
    Ai::Float64                    # what is this for????
    b_friction::Float64
    c_friction::Float64
    Rb::Float64
    lim_dz::Int
    mtq_max::Vector{Float64}
    rw_max_torque::Float64
    p_500::Float64
    Cm::Vector{Float64}
    known_rm::Vector{Float64}
end

function create_satellite_parameters()
    orbits = 8
    # m = 4
    # lx = 0.1
    # ly = 0.1
    # lz = 0.3405

    N_coils = [400, 400, 800]
    A_coils = [0.0022, 0.0022, 0.001566]
    R_coils = [110, 110, 31.32]
    mt = [5, 5, 1]
    W_up_limit = [0.2273, 0.2273, 0.8017]
    I_up_limit = [0.0455, 0.0455, 0.1566]

    Km = 20e-5
    Rb = 10
    Ai = 1

    Kv = 1 / Km
    b_friction = 9.5e-9
    c_friction = 1.9e-7
    Jw = 1.9e-6
    A = 0.12    
    lim_dz = 300

    p_500 = 1.80e-12

    # Ix = (m / 12) * (ly^2 + lz^2)
    # Iy = (m / 12) * (lx^2 + lz^2)
    # Iz = (m / 12) * (lx^2 + ly^2)

    PAI = [-1 0.01 -0.04; 0.01 1 0.00; -0.04 0.00 1.00]
    PMI = diagm(0 => [0.04127073921, 0.041018997570, 0.00690030456])
    Cm = [0.00415, 0.00116, 0.0016]

    for j = 1:3
        PAI[:, j] = PAI[:, j] / norm(PAI[:, j])
    end
    I = PAI * PMI * PAI'
    I_inv = inv(I)                  # Here I just used the inv() function, instead of this formula I_inv = eye(3,3)/I

    Re = 6371.2e3
    Rs = 500e3
    Radius = Re + Rs
    G = 6.67428e-11
    M = 5.972e24
    w_o = sqrt(G * M / Radius^3)
    v_satellite = sqrt(G * M / Radius)
    w_o_io = [0, w_o, 0]

    known_rm = [0.048, 0.051, 0.047]

    orbitPeriod = (2 * pi) / w_o
    n = orbitPeriod * orbits
    mtq_max = [0.2, 0.2, 0.2]
    rw_max_torque = 1e-4

    return SatelliteParameters(Radius,
                               w_o,
                               v_satellite, 
                               n, 
                               orbitPeriod, 
                               I, 
                               I_inv,
                               w_o_io, 
                               N_coils, 
                               A_coils, 
                               R_coils, 
                               mt, 
                               I_up_limit, 
                               W_up_limit, 
                               Jw, 
                               A, 
                               Km, 
                               Kv, 
                               Ai, 
                               b_friction, 
                               c_friction, 
                               Rb, 
                               lim_dz, 
                               mtq_max, 
                               rw_max_torque, 
                               p_500, 
                               Cm, 
                               known_rm)
end
