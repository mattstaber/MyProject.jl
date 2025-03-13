using LinearAlgebra
using DifferentialEquations
using CairoMakie
using Statistics

# Include all Julia files
include("aux_dyn.jl")
include("aux_dynE_2D.jl")
include("aux_dynE.jl")
include("cal_airdens.jl")
include("cal_BR_prdt.jl")
include("cal_BR.jl")
include("cal_delpsi.jl")
include("cal_drdc.jl")
include("cal_J.jl")
include("cal_pci2pcr.jl")
include("cal_sphdis.jl")
include("cal_z_3d.jl")
include("cal_z.jl")
include("myplot.jl")
include("NPCG_2024.jl")
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


    # Normalization
    DU = rp
    VU = sqrt(mu / DU)
    TU = DU / VU
    rp = rp / DU
    mu = mu * DU^-3 * TU^2
    w = w * TU
    Vf = Vf / VU
    rf = rf / DU
    V0 = V0 / VU
    r0 = r0 / DU

    # Initial state vector
    X0 = [r0, theta0, phi0, V0, gamma0, psi0]

    # Initial and final energy conditions
    e0 = mu / r0 - (V0^2) / 2
    ef = mu / rf - (Vf^2) / 2

    # Initial bank reversal value
    delpsi0 = cal_delpsi(X0, thetatgt, phitgt)
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
        "rp" => rp,
        "mu" => mu,
        "w" => w,
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
        "ef" => ef,
        "rf" => rf,
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
    tc = 0  # current time
    ec = e0  # current energy
    hc = r0 - rp  # current altitude
    Vc = V0  # current velocity
    sigcmdprv = 0  # previous step bank angle command
    sigcmdprv2 = 0  # two-step previous bank angle command

    X = Vector{Vector{Float64}}()
    T = Vector{Float64}()
    SIG = Vector{Float64}()

    push!(X, X0)      # Xn is already a Float64 array
    push!(T, tc)
    push!(SIG, 0.0)  # bankcmd should be a Float64 value

    testdata = Vector{Vector{}}()

    # ------------------- software-in-the-loop Simulation ------------------- #
    tic()
    while ec < ef
        bankcmd, sigs, BRprev, sigcmdprv, sigcmdprv2, testdata = npcg(tc, ec, TU, GAT, DU, siglmt, tstp, sigdlmt, sigddlmt, sigcmdprv, sigcmdprv2, BRprev, sigs, BAL, BRL, X0, thetatgt, phitgt, Case_Number, auxdata, testdata)

        # Time span for ODE solver
        tspan = (tc, tc + tstp / TU)

        # parameters
        p = [bankcmd, auxdata]

        # Define the ODEProblem
        prob = ODEProblem(aux_dyn!, X0, tspan, p)

        # Solve the ODE
        sol = solve(prob, BS5(), saveat=0.1, reltol=1e-6, abstol=1e-9)

        # Extract time points and solution
        tt = sol.t  # Time points
        XX = sol.u  # Solution at those time points

        Xn = XX[end]
        tc = tt[end]

        # Update for the next step
        hc = Xn[1] - rp
        Vc = Xn[4]
        gamc = Xn[5]
        X0 = Xn
        ec = mu / Xn[1] - (Xn[4]^2) / 2

        # Save data
        push!(X, Xn)
        push!(T, tc)
        push!(SIG, bankcmd)

        # Check guidance termination
        Rgo = DU * cal_sphdis(X0[2], X0[3], thetatgt, phitgt)
        if Rgo < 100  # meters
            println("Rgo < 100 m")
            break
        end
    end
    elpdtime = toc();

    # --------------------------- post-processing --------------------------- #
    r = [x[1] for x in X]
    theta = [x[2] for x in X]
    phi = [x[3] for x in X]
    V = [x[4] for x in X]
    gamma = [x[5] for x in X]
    psi = [x[6] for x in X]

    # Energy-like variable
    e = mu ./ r - 0.5 * V .^ 2

    # Time conversion
    time = T .* TU

    # Bank angle
    sigcmd = SIG
    sigmadata = hcat(time, sigcmd)

    # Downrange and crossrange initialization
    nn = length(time)
    DR = zeros(nn)
    CR = zeros(nn)
    delpsi = zeros(nn)

    for ii in 1:nn
        dr, cr = cal_drdc(theta0, phi0, thetatgt, phitgt, theta[ii], phi[ii])
        DR[ii] = dr
        CR[ii] = cr
        delpsi[ii] = cal_delpsi(X[ii], thetatgt, phitgt)
    end

    DR0 = cal_sphdis(theta0, phi0, thetatgt, phitgt)  # Mission range-to-go
    DRtogo = (DR0 .- DR) .* DU  # Downrange-to-go
    Rgof = cal_sphdis(theta[end], phi[end], thetatgt, phitgt) .* DU .* 1e-3

    # Path constraints
    Vreal = V .* VU  # m/s
    hreal = (r .- rp) .* DU  # m
    rho = cal_airdens(hreal, pn)  # kg/m^3
    q = 0.5 .* rho .* Vreal .^ 2  # Pa
    D = q .* CD .* S ./ m  # m/s^2
    L = q .* CL .* S ./ m  # m/s^2
    A = sqrt.(L .^ 2 .+ D .^ 2) ./ 9.81  # g
    qdot = kq .* (rho .^ kqN) .* (Vreal .^ kqM)  # kW/m^2

    # Unit conversion
    h = (r .- rp) .* DU .* 1e-3
    V = V .* VU .* 1e-3
    theta = theta .* r2d
    phi = phi .* r2d
    gamma = gamma .* r2d
    psi = psi .* r2d
    DR = DR .* DU .* 1e-3
    CR = CR .* DU .* 1e-3
    DRtogo = DRtogo .* 1e-3

    # Final targeting error and simulation time
    println("\nSimulation Time: ", round(elpdtime, digits=3), " sec")

    f = plot_edl_data(sigmadata, V, h, DR, theta, phi, time, gamma, psi, CR, A, q, qdot, r2d)
    display(f)

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
Case_Number = 5
BAP = 1

run_simulation(Case_Number,BAP)