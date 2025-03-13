# ----------------------------------------------------------------------- #
# problem formulation for the atmospheric entry phase of the
# human Mars EDL mission (Mid lift-to-drag ratio rigid vehicle (MRV))
# Main reference: Sagliano et al., 2024-1171 (AIAA SCITECH 2024 Forum)
#
#    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
#    Prof. Dae Young Lee** and Prof. Bong Wie***
#    Iowa State University, Ames, IA 50011
#    *leeyr111@gmail.com
#    **daylee@iastate.edu
#    ***bongwie@iastate.edu
# ----------------------------------------------------------------------- #
function ex_Mars_manned_new(BAP)
    d2r = pi / 180
    r2d = 180 / pi

    # Constants
    rp = 3396.19e3       # Mars radius, m
    mu = 4.2828e13       # Gravitational parameter, m^3/s^2
    w = 7.0882e-5        # Spin rate, rad/s

    # Planet selection (1: Earth, 2: Mars)
    pn = 2

    # Vehicle parameters (L/D = 0.54, B = 379 kg/m^2)
    S = 160              # Cross-sectional area, m^2
    m = 58800            # Mass, kg
    CL = 0.5236         # Lift coefficient
    CD = 0.9696         # Drag coefficient

    # Entry interface in planet-centered-rotating frame
    r0 = rp + 125e3      # Radius, m
    theta0 = -176.40167 * d2r  # Longitude, rad
    phi0 = -21.3 * d2r   # Latitude, rad
    V0 = 4700           # Velocity, m/s
    gamma0 = -10 * d2r  # Flight path angle, rad
    psi0 = -2.8758 * d2r  # Heading angle, rad

    # Final conditions
    hf = 2480
    rf = hf + rp
    Vf = 450
    thetatgt = -175.8 * d2r
    phitgt = 0.276 * d2r

    # Path constraints
    Amax = 4           # g-load, g
    qmax = 13          # Dynamic pressure, kPa
    Qdotmax = 500      # Heating rate, kW/m^2

    # Heating rate constants
    kq = 5.3697e-8
    kqN = 0.5
    kqM = 3.15

    # Bank angle constraints
    BAL = 1
    siglmt = 180 * d2r   # Bank magnitude limit
    sigdlmt = 20 * d2r   # Bank rate limit
    sigddlmt = 20 * d2r  # Bank acceleration limit

    # Guidance activation time
    GAT = 170

    # Bank reversal logic
    BRL = 2
    dlpsT = [6, 3] .* d2r  # Convert to radians
    KBR = 7

    # Bank angle parameterization options
    sig0 = 70 * d2r  # Initial guess
    sigf = 0
    KEF = 0
    KLF = 0
    
    if BAP == 1
        sigf = 40 * d2r
    elseif BAP == 2
        sigf = 10 * d2r
        KEF = 1.1
    elseif BAP == 3
        KEF = 0.6
    elseif BAP == 4
        KLF = 1.3
    else
        error("Invalid BAP value")
    end
    sigs = sig0

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
