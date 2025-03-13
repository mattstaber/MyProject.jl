function ex_Mars_robotic(BAP)
    # Constants
    rp = 3397e3    # Mars radius, m
    mu = 4.2828e13 # Gravitational parameter, m^3/s^2
    w  = 7.0882e-5 # Spin rate, rad/s
    d2r = π / 180
    r2d = 180 / π

    # Planet selection: 1 = Earth, 2 = Mars
    pn = 2

    # Vehicle parameters (L/D = 0.24, B = 146 kg/m²)
    S = 15.588  # Cross-sectional area, m²
    m = 3300    # Mass, kg
    CL = 0.348  # Lift coefficient
    CD = 1.45   # Drag coefficient

    # Entry interface in planet-centered-inertial frame (Dutta and Braun 2014)
    r0     = rp + 125e3   # Radius, m
    theta0 = 126.72 * d2r # Longitude, rad
    phi0   = -3.9186 * d2r # Latitude, rad
    V0     = 6083.3       # Velocity, m/s
    gamma0 = -15.4892 * d2r # Flight path angle, rad
    psi0   = 93.2065 * d2r  # Heading angle, rad

    # Coordinate transformation: Inertial frame to rotating frame
    r0, theta0, phi0, V0, gamma0, psi0 = 
        cal_pci2pcr(r0, theta0, phi0, V0, gamma0, psi0, w)

    # Final condition
    hf       = 12000
    rf       = hf + rp
    Vf       = 406
    thetatgt = 137.431 * d2r
    phitgt   = -4.552 * d2r

    # Path constraints limits
    Amax    = 15    # G-load, g
    qmax    = 21    # Dynamic pressure, kPa
    Qdotmax = 1000  # Heating rate, kW/m²

    # Heating rate constants
    kq  = 5.3697e-8
    kqN = 0.5
    kqM = 3.15

    # Bank angle constraints
    BAL     = 1
    siglmt  = 180 * d2r  # Bank magnitude limit
    sigdlmt = 15 * d2r   # Bank rate limit
    sigddlmt = 15 * d2r  # Bank acceleration limit

    # Guidance activation time
    GAT = 60

    # Bank reversal logic
    BRL   = 2
    dlpsT = [4, 2] * d2r
    KBR   = 8

    # Bank angle parameterization options
    auxdata = Dict()
    KEF = 0
    KLF = 0
    if BAP == 1
        sig0 = 70 * d2r  # Initial guess
        sigf = 20 * d2r  # Design parameter
        sigs = sig0
        auxdata[:sigf] = sigf

    elseif BAP == 2
        sig0 = 70 * d2r
        sigf = 10 * d2r
        KEF  = 1.2
        sigs = sig0
        auxdata[:sigf] = sigf
        auxdata[:KEF] = KEF

    elseif BAP == 3
        sig0 = 70 * d2r
        KEF  = 1.0
        sigs = sig0
        auxdata[:KEF] = KEF

    elseif BAP == 4
        sig0 = 70 * d2r
        KLF  = 1.4
        sigs = sig0
        auxdata[:KLF] = KLF
    else
        error("Invalid BAP value. Choose 1, 2, 3, or 4.")
    end

    # Return all relevant data as a dictionary
    return Dict(
        :rp => rp, :mu => mu, :w => w, :d2r => d2r, :r2d => r2d, :pn => pn,
        :S => S, :m => m, :CL => CL, :CD => CD,
        :r0 => r0, :theta0 => theta0, :phi0 => phi0, :V0 => V0,
        :gamma0 => gamma0, :psi0 => psi0,
        :hf => hf, :rf => rf, :Vf => Vf, :thetatgt => thetatgt, :phitgt => phitgt,
        :Amax => Amax, :qmax => qmax, :Qdotmax => Qdotmax,
        :BTU2kW => 0, :kq => kq, :kqN => kqN, :kqM => kqM,
        :BAL => BAL, :siglmt => siglmt, :sigdlmt => sigdlmt, :sigddlmt => sigddlmt,
        :GAT => GAT, :BRL => BRL, :dlpsT => dlpsT, :KBR => KBR,
        :sig0 => sig0, :sigf => sigf, :sigs => sigs, :KEF => KEF, :KLF => KLF
    )
end