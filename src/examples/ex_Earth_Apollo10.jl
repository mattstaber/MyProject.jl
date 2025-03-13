# ----------------------------------------------------------------------- #
# problem formulation for the atmospheric entry phase of the
# human Earth EDL mission, Apollo 10 command module
# Main reference: Szelc 1969, NASA TM-X-69555
#
#    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
#    Prof. Dae Young Lee** and Prof. Bong Wie***
#    Iowa State University, Ames, IA 50011
#    *leeyr111@gmail.com
#    **daylee@iastate.edu
#    ***bongwie@iastate.edu
# ----------------------------------------------------------------------- #

function ex_Earth_Apollo10(BAP)
    d2r = pi / 180
    r2d = 180 / pi
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # Problem formulation for the atmospheric entry phase of the
    # human Earth EDL mission, Apollo 10 command module
    # Main reference: Szelc 1969, NASA TM-X-69555

    # Constants
    rp = 6378.137e3      # Earth radius, m
    mu = 398600.435507e9  # Gravitational parameter, m^3/s^2
    w = 2 * π / 86400     # Spin rate, rad/s

    # Planet / 1: Earth, 2: Mars
    pn = 1

    # Vehicle parameters
    S = 12.017     # Cross-sectional area, m^2
    m = 5498.22    # Mass, kg
    CL = 0.40815   # Lift coefficient
    CD = 1.2569    # Drag coefficient

    # Entry interface in planet-centered-inertial frame (Szelc 1969)
    r0 = 6498.270e3             # Geocentric radius, m
    theta0 = 174.24384 * d2r    # Longitude, deg
    phi0 = -23.51457 * d2r      # Geocentric latitude, m
    V0 = 11.06715e3             # Inertial velocity, m/s
    gamma0 = -6.6198381 * d2r   # Inertial flight path angle
    psi0 = 71.9317 * d2r        # Inertial heading angle

    # Coordinate transformation, inertial frame to rotating frame
    r0, theta0, phi0, V0, gamma0, psi0 = cal_pci2pcr(r0, theta0, phi0, V0, gamma0, psi0, w)

    # Final condition
    hf = 5683
    rf = hf + rp
    Vf = 119.67
    thetatgt = 194.8501 * d2r
    phitgt = -15.4121 * d2r

    # Path constraints (Szelc 1969, Pavlosky and Leger 1974)
    Amax = 10       # g-load, g
    qmax = 30       # Dynamic pressure, kPa
    Qdotmax = 6000  # Heating rate, kW/m^2

    # Heating rate constant (Young and Smith 1967)
    BTU2kW = 11.356539  # 1 BTU/s.ft2 to kW/m^2
    kq = BTU2kW * 20 * 1e-9
    kqN = 0.5
    kqM = 3

    # Bank angle constraints
    BAL = 1
    siglmt = 180 * d2r   # Bank magnitude limit
    sigdlmt = 15 * d2r   # Bank rate limit
    sigddlmt = 15 * d2r  # Bank acceleration limit

    # Guidance activation time
    GAT = 90

    # Bank reversal logic
    BRL = 1
    dlpsT = [5, 7] .* d2r  # Convert to radians
    KBR = 6

    # Bank angle parameterization options
    sig0 = 70 * d2r # initial guess of the parameter to be found
    KEF = 0
    KLF = 0
    sigf = 0
    if BAP == 1
        sigf = 30 * d2r # design parameter
    elseif BAP == 2
        sigf = 10 * d2r # design parameter
        KEF = 0.9       # design parameter
    elseif BAP == 3
        KEF = 1.3       # design parameter
    elseif BAP == 4
        KLF = 1.3       # design parameter
    else
        error("Invalid BAP value")
    end
    sigs = sig0

    # Return all variables as a Dict (or NamedTuple for structure)
    return Dict(
        :rp => rp, :mu => mu, :w => w, :d2r => d2r, :r2d => r2d, :pn => pn,
        :S => S, :m => m, :CL => CL, :CD => CD,
        :r0 => r0, :theta0 => theta0, :phi0 => phi0, :V0 => V0,
        :gamma0 => gamma0, :psi0 => psi0,
        :hf => hf, :rf => rf, :Vf => Vf, :thetatgt => thetatgt, :phitgt => phitgt,
        :Amax => Amax, :qmax => qmax, :Qdotmax => Qdotmax,
        :BTU2kW => BTU2kW, :kq => kq, :kqN => kqN, :kqM => kqM,
        :BAL => BAL, :siglmt => siglmt, :sigdlmt => sigdlmt, :sigddlmt => sigddlmt,
        :GAT => GAT, :BRL => BRL, :dlpsT => dlpsT, :KBR => KBR,
        :sig0 => sig0, :sigf => sigf, :sigs => sigs, :KEF => KEF, :KLF => KLF
    )
end