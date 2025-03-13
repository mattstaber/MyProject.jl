
# ----------------------------------------------------------------------- #
# problem formulation for the atmospheric entry phase of the
# robotic Earth EDL mission, CAV-H
# Main reference: Lu, AIAA GNC conference 2013-4648
#
#    Coded by Dr. Youngro Lee* for his PhD work under the supervision of
#    Prof. Dae Young Lee** and Prof. Bong Wie***
#    Iowa State University, Ames, IA 50011
#    *leeyr111@gmail.com
#    **daylee@iastate.edu
#    ***bongwie@iastate.edu
# ----------------------------------------------------------------------- #
function ex_Earth_CAVH(BAP)
    # constants
    rp = 6378.137e3      # Earth radius, m
    mu = 398600.435507e9 # gravitational parameter, m^3/s^2
    w = 2 * pi / 86400      # spin rate, rad/s
    d2r = pi / 180
    r2d = 180 / pi

    # planet / 1: Earth, 2: Mars
    pn = 1

    # vehicle parameters (Phillips 2003)
    S = 0.4839   # cross-sectional area, m^2
    m = 907      # mass, kg
    CL = 0.42324 # lift cofficient, mean value
    CD = 0.1292  # drag cofficient, mean value

    # entry interface in planet-centered-rotating frame (Lu 2013)
    # no need to tranform the coordinate
    r0 = rp + 122e3
    theta0 = -72.2744 * d2r
    phi0 = 39.1814 * d2r
    V0 = 7400
    gamma0 = -1 * d2r
    psi0 = 38.5668 * d2r

    # final condition
    hf = 30e3
    rf = hf + rp
    Vf = 2000
    thetatgt = 64.1166 * d2r
    phitgt = 31.4927 * d2r

    # path constraints
    Amax = 2.5  # g-load, g
    qmax = 100  # dynamic pressure, kPa
    Qdotmax = 8000 # heating rate, kW/m^2

    # heating rate constant
    kq = 9.4369e-8
    kqN = 0.5
    kqM = 3.15

    # bank angle constraints
    BAL = 1
    siglmt = 90 * d2r   # bank magnitude limit
    sigdlmt = 15 * d2r  # bank rate limit
    sigddlmt = 15 * d2r  # bank acc limit

    # guidance activation time
    GAT = 300

    # bank reversal logic
    BRL = 1
    dlpsT = [15 5] * d2r
    KBR = 6

    # bank angle parameterization options
    BTU2kW = 0
    KEF = 0
    KLF = 0
    sigf = 0
    if BAP == 1

        sig0 = 70 * d2r # initial guess of the parameter to be found
        sigf = 60 * d2r # design parameter

        sigs = sig0

    elseif BAP == 2

        sig0 = 70 * d2r # initial guess of the parameter to be found
        sigf = 60 * d2r # design parameter
        KEF = 0.9    # design parameter

        sigs = sig0


    elseif BAP == 3

        sig0 = 70 * d2r # initial guess of the parameter to be found

        KEF = 0.6 # design parameter

        sigs = sig0

    elseif BAP == 4

        sig0 = 70 * d2r # initial guess of the parameter to be found

        KLF = 0.5 # design parameter

        sigs = sig0

    end

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