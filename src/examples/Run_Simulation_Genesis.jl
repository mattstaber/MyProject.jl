using LinearAlgebra
using DifferentialEquations
using CairoMakie
using Genesis

# Include all Julia files
# include("aux_dyn.jl")
# include("aux_dynE_2D.jl")
# include("aux_dynE.jl")
# include("cal_airdens.jl")
# include("cal_BR_prdt.jl")
# include("cal_BR.jl")
# include("cal_delpsi.jl")
# include("cal_drdc.jl")
# include("cal_J.jl")
# include("cal_pci2pcr.jl")
# include("cal_sphdis.jl")
# include("cal_z_3d.jl")
# include("cal_z.jl")
# include("/src/functions/NPCG.jl")
include("ex_Earth_Apollo10.jl")
include("ex_Earth_CAVH.jl")
include("ex_Mars_manned.jl")
include("ex_Mars_manned_new.jl")
include("ex_Mars_robotic.jl")

function tic()
    global start_time = Base.time()
end

function toc()
    return Base.time() - start_time
end

struct MyAtmosphere{S<:Vehicle} <: Atmosphere
    state::S
end

function Genesis.density(atmos::MyAtmosphere)
    return cal_dens(planetodetic_altitude(atmos.state))
end

function cal_dens(h)
    return 1.225 * exp(-0.14(h / 1000))
end

struct MarsAtmosphere{S<:Vehicle} <: Atmosphere
    state::S
end

function Genesis.density(atmos::MarsAtmosphere)
    return cal_dens_Mars(planetodetic_altitude(atmos.state))
end

function cal_dens_Mars(h)
    # atmosphere (AAS 18-485)
    beta1 = 559.351005946503
    beta2 = 188.95110711075

    T = (1.4e-13) .* (h .^ 3) .- (8.85e-9) .* (h .^ 2) .- (1.245e-3) .* h .+ 205.3645
    rho0 = beta1 ./ (beta2 .* T)
    beta = -0.000105 .* h
    rho = rho0 .* exp.(beta)  # Use broadcasting (.)
    return rho
    #return 1.225 * exp(-0.14(h / 1000))
end

function run_simulation(Case_Number, BAP)
    # Constants
    d2r = π / 180.0  # Degrees to radians conversion factor
    r2d = 180 / π

    # ----------------------------------------------------------------------- #
    # Entry example cases
    if Case_Number == 1
        result = ex_Mars_robotic(BAP)
    elseif Case_Number == 2
        result = ex_Mars_manned(BAP)
    elseif Case_Number == 3
        result = ex_Mars_manned_new(BAP)
    elseif Case_Number == 4
        result = ex_Earth_Apollo10(BAP)
    elseif Case_Number == 5
        result = ex_Earth_CAVH(BAP)
    else
        println("Invalid case")
    end

    # Unpack the dictionary into variables
    rp = result[:rp]
    mu = result[:mu]
    w = result[:w]
    d2r = result[:d2r]
    r2d = result[:r2d]
    pn = result[:pn]

    S = result[:S]
    m = result[:m]
    CL = result[:CL]
    CD = result[:CD]

    r0 = result[:r0]
    theta0 = result[:theta0]
    phi0 = result[:phi0]
    V0 = result[:V0]
    gamma0 = result[:gamma0]
    psi0 = result[:psi0]

    hf = result[:hf]
    rf = result[:rf]
    Vf = result[:Vf]
    thetatgt = result[:thetatgt]
    phitgt = result[:phitgt]

    Amax = result[:Amax]
    qmax = result[:qmax]
    Qdotmax = result[:Qdotmax]

    BTU2kW = result[:BTU2kW]
    kq = result[:kq]
    kqN = result[:kqN]
    kqM = result[:kqM]

    BAL = result[:BAL]
    siglmt = result[:siglmt]
    sigdlmt = result[:sigdlmt]
    sigddlmt = result[:sigddlmt]

    GAT = result[:GAT]
    BRL = result[:BRL]
    dlpsT = result[:dlpsT]
    KBR = result[:KBR]

    sig0 = result[:sig0]
    sigf = result[:sigf]
    sigs = result[:sigs]
    KEF = result[:KEF]
    KLF = result[:KLF]
    # ----------------------------------------------------------------------- #

    time = Genesis.Time()
    integ = Genesis.DP5(time, max_dt=1)

    executive = Executive(integ,
        stop_time=5140)

    planet = Planet(PointMassGravity(gravitational_parameter=mu),
        PlanetOrientation(rotation_rate=w),
        Ellipsoid(equatorial_radius=rp, flattening=0),
        time)

    vehicle = Vehicle(planet)

    add!(executive, vehicle)

    set_mass_properties!(vehicle,
        mass=m)

    configure!(vehicle,
        ComputedAcceleration(vehicle))

    set_position!(vehicle,
        radius=r0,
        longitude=theta0,
        declination=phi0)

    set_velocity!(vehicle,
        planet_relative_velocity=V0,
        planetodetic_planet_relative_flight_path_angle=gamma0,
        planetodetic_planet_relative_azimuth=psi0)

    configure!(vehicle,
        PrescribedAttitude(vehicle,
            angle_of_attack=0.0 / 180 * pi,
            sideslip_angle=0.0,
            bank_angle=0 * pi / 180))

    if Case_Number == 1
        configure!(vehicle,
        MarsAtmosphere(vehicle))
    elseif Case_Number == 2
        configure!(vehicle,
        MarsAtmosphere(vehicle))
    elseif Case_Number == 3
        configure!(vehicle,
        MarsAtmosphere(vehicle))
    elseif Case_Number == 4
        configure!(vehicle,
        MyAtmosphere(vehicle))
    elseif Case_Number == 5
        configure!(vehicle,
        MyAtmosphere(vehicle))
    else
        println("Invalid case")
    end

    configure!(vehicle,
        MyAtmosphere(vehicle))

    configure!(vehicle,
        NonsymmetricAerodynamics(vehicle,
            reference_area=S,
            cd=CD,
            cl=CL,
            cs=0))

    # Normalization
    DU = rp
    VU = sqrt(mu / DU)
    TU = DU / VU
    rp_norm = rp / DU
    mu_norm = mu * DU^(-3) * TU^2
    w_norm = w * TU
    Vf_norm = Vf / VU
    rf_norm = rf / DU
    V0_norm = V0 / VU
    r0_norm = r0 / DU

    # Initial state vector
    X0_norm = [r0_norm, theta0, phi0, V0_norm, gamma0, psi0]

    # Initial and final energy conditions
    e0_norm = mu_norm / r0_norm - (V0_norm^2) / 2
    ef_norm = mu_norm / rf_norm - (Vf_norm^2) / 2

    # Initial bank reversal value
    delpsi0 = cal_delpsi(X0_norm, thetatgt, phitgt)

    if sign(delpsi0) > 0
        BRprev = -1
    else
        BRprev = 1
    end

    # Increment for Newton-Raphson method
    dsig = 0.01 * d2r

    # Propagation time step and guidance frequency
    tstp = 1  # s

    auxdata = Dict(
        "pn" => pn,
        "rp" => rp_norm,
        "mu" => mu_norm,
        "w" => w_norm,
        "S" => S,
        "m" => m,
        "CL" => CL,
        "CD" => CD,
        "DU" => DU,
        "VU" => VU,
        "TU" => TU,
        "theta0" => theta0,
        "phi0" => phi0,
        "thetatgt" => thetatgt,
        "phitgt" => phitgt,
        "ef" => ef_norm,
        "rf" => rf_norm,
        "BAP" => BAP,
        "KBR" => KBR,
        "dlpsT" => dlpsT,
        "opt" => 0,
        "dsig" => dsig,
        "sigf" => sigf,
        "KEF" => KEF,
        "KLF" => KLF
    )

    # Initialization
    tc = 0
    ec_norm = e0_norm  # current energy
    sigcmdprv = 0  # previous step bank angle command
    sigcmdprv2 = 0  # two-step previous bank angle command

    # ------------------- software-in-the-loop Simulation ------------------- #
    configure!(executive, FlightSoftware([FlightSoftwareRateGroup(dt=1 // 1) do
        tc = dynamic_time(time)

        X0 = [planetodetic_altitude(vehicle) + rp, mod(longitude(vehicle) + 2pi, 2pi), planetodetic_latitude(vehicle), norm(planet_relative_velocity(vehicle, frame=PCPF(planet))), planetodetic_planet_relative_flight_path_angle(vehicle), planetodetic_planet_relative_azimuth(vehicle)]
        X0_norm = [X0[1] / DU, X0[2], X0[3], X0[4] / VU, X0[5], X0[6]]
        ec_norm = (mu_norm) / X0_norm[1] - (X0_norm[4]^2) / 2

        bankcmd, sigs, BRprev, sigcmdprv, sigcmdprv2 = npcg(tc / TU, ec_norm, TU, GAT, DU, siglmt, tstp, sigdlmt, sigddlmt, sigcmdprv, sigcmdprv2, BRprev, sigs, BAL, BRL, X0_norm, thetatgt, phitgt, Case_Number, auxdata)



        configure!(vehicle,
            PrescribedAttitude(vehicle,
                angle_of_attack=0.0 / 180 * pi,
                sideslip_angle=0.0,
                bank_angle=bankcmd))
    end]))

    add!(executive,
        Event(condition=() -> ec_norm * 1.0001 > ef_norm,
            action=() -> stop!(executive)))

    time_history = TimeHistoryLoggingGroup(
        "t" => () -> dynamic_time(time),
        "h" => () -> planetodetic_altitude(vehicle),
        "vi" => () -> inertial_velocity(vehicle, frame=PCPF(planet)),
        "vr" => () -> planet_relative_velocity(vehicle, frame=PCPF(planet)),
        "ha" => () -> radius_of_apoapsis(vehicle) - equatorial_radius(planet),
        "ε" => () -> specific_energy(vehicle),
        "lat" => () -> planetodetic_latitude(vehicle) * r2d,
        "lon" => () -> mod(longitude(vehicle) + 2pi, 2pi) * r2d,
        "gam" => () -> planetodetic_planet_relative_flight_path_angle(vehicle) * r2d,
        "psi" => () -> planetodetic_planet_relative_azimuth(vehicle),
        "sig" => () -> bank_angle(vehicle),
        "dp" => () -> dynamic_pressure(vehicle) / 1000,
        "acc" => () -> norm(acceleration(vehicle, frame=PCI(planet))) / 9.81,
        "drag" => () -> dynamic_pressure(vehicle) * S * CD,
        "lift" => () -> dynamic_pressure(vehicle) * S * CL,
        "gload" => () -> sqrt((dynamic_pressure(vehicle) / 364.02) .^ 2 + (dynamic_pressure(vehicle) / 364.02 * CL) .^ 2) / 9.81
    )

    add!(executive, time_history)

    tic()
    run!(executive)
    elpdtime = toc()

    # Final targeting error and simulation time
    println("\nSimulation Time: ", round(elpdtime, digits=3), " sec")

    outputs = data(time_history)

    f = Figure(size=(600, 600))
    ax = Axis(f[1, 1],
        xlabel="Time [s]",
        ylabel="Bank angle [deg]")
    lines!(ax,
        outputs["t"],
        outputs["sig"] / pi * 180)

    ax = Axis(f[1, 2],
        xlabel="Planet-Relative Velocity [km/s]",
        ylabel="Altitude [km]")
    lines!(ax,
        norm.(outputs["vr"]) ./ 1000,
        outputs["h"] ./ 1000)

    ax = Axis(f[2, 1],
        xlabel="Time [s]",
        ylabel="Flight path angle [deg]")
    lines!(ax,
        outputs["t"],
        outputs["gam"])

    ax = Axis(f[2, 2],
        xlabel="Longitude [deg]",
        ylabel="Latitude [deg]")
    lines!(ax,
        outputs["lon"],
        outputs["lat"])

    ax = Axis(f[3, 1],
        xlabel="Time [s]",
        ylabel="g-load [g]")
    lines!(ax,
        outputs["t"],
        outputs["gload"])
    ax = Axis(f[3, 2],
        xlabel="Time [s]",
        ylabel="Dynamic pressure [kPa]")
    lines!(ax,
        outputs["t"],
        outputs["dp"])
    f
end



# ----------------------------------------------------------------------- #
#  select the entry example (Case_Number) out of five cases
#  select the bank angle parameterization (BAP) out of seven options
# ----------------------------------------------------------------------- #
# Case_Number = atmospheric entry examples
#     1: Mars robotic mission (Mars Science Laboratory model)
#     2: Mars manned mission model in (Jiang et al., 2019)
#     3: Mars manned mission model updated in (Signialo et al, 2024)
#     4: Apollo 10 (Szelc 1969)
#     5: CAV-H (Lu 2013)
#
# BAP = bank angle parameterization
#     1: linear function            (Lu 2014)
#     2: exponential function       (Liang and Zhu 2021)
#     3: exponential function       (Youngro Lee)
#     4: logistic function          (Youngro Lee)
#
# Note!!! each BAP option has guidance parameters to be adjusted, which
#         can be found in "ex_" scripts.
#
# ----------------------------------------------------------------------- #
#  GAT, BRL, and BAL should be wisely selected for a good guidance
#  performance. A good parameter set for each combination are already given
#  in "ex_***", and they can be adjusted according to the user's desire.
# ----------------------------------------------------------------------- #
#
# GAT = guidance activation time in seconds since entry interface (EI)
#
# BRL = bank reversal logic
#     1: conventional logic by Tu, Munir, Mease, and Bayard (JGCD 2000)
#           dlpsT = heading angle error thresholds
#     2: predictive logic by K.M. Smith (AAS 2016)
#           KBR = damping ratio of the predictive lateral guidance
#
# BAL = bank angle constraints
#     1: simple magnitude contraint
#           siglmt   = magnitude limit
#     2: rate and acceleration constraints
#           sigdlmt  = rate limit
#           sigddlmt = acceleration limit
#
# ----------------------------------------------------------------------- #
Case_Number = 3
BAP = 1

run_simulation(Case_Number, BAP)