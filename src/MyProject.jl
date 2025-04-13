module MyProject

using Genesis
using StaticArrays
using LinearAlgebra
using DifferentialEquations

mutable struct GuidanceState
    sigs::Float64
    BR::Int
    sigcmdprv::Float64
    sigcmdprv2::Float64
    first_call::Bool
end

function _normalize(vehicle, mu, rf, Vf)
    DU = equatorial_radius(planet(vehicle))
    VU = sqrt(mu / DU)
    TU = DU / VU
    mu_norm = mu * DU^(-3) * TU^2
    rf_norm = rf / DU
    Vf_norm = Vf / VU

    return DU, VU, TU, mu_norm, rf_norm, Vf_norm
end

function npcg_step!(
    state::GuidanceState,
    vehicle::Vehicle,
    theta0::Float64,
    phi0::Float64,
    thetatgt::Float64,
    phitgt::Float64,
    rf::Float64,
    Vf::Float64,
    GAT::Int64,
    bank_constraints::NamedTuple,
    KBR::Int64,
    guidance_logic::NamedTuple,
)::Float64

    p = planet(vehicle)
    mu = gravitational_parameter(p)
    w = angular_velocity(p; frame=PCPF(p))[3]

    DU, VU, TU, mu_norm, rf_norm, Vf_norm = _normalize(vehicle, mu, rf, Vf)
    tc = dynamic_time(Genesis.time(p)) / TU

    X0 = [
        planetodetic_altitude(vehicle) + equatorial_radius(p),
        mod(longitude(vehicle) + 2π, 2π),
        planetodetic_latitude(vehicle),
        norm(planet_relative_velocity(vehicle, frame=PCPF(p))),
        planetodetic_planet_relative_flight_path_angle(vehicle),
        planetodetic_planet_relative_azimuth(vehicle)
    ]

    X0_norm = [X0[1] / DU, X0[2], X0[3], X0[4] / VU, X0[5], X0[6]]
    ec_norm = mu_norm / X0_norm[1] - (X0_norm[4]^2) / 2
    ef_norm = mu_norm / rf_norm - (Vf_norm^2) / 2

    auxdata = Dict(
        "rp" => equatorial_radius(p) / DU,
        "mu" => mu_norm,
        "w" => w * TU,
        "S" => Genesis.reference_area(aerodynamics(vehicle)),
        "m" => Genesis.element_mass(vehicle),
        "CL" => Genesis.force_coefficients(aerodynamics(vehicle)).cl.undispersed,
        "CD" => Genesis.force_coefficients(aerodynamics(vehicle)).cd.undispersed,
        "DU" => DU,
        "VU" => VU,
        "TU" => TU,
        "theta0" => theta0,
        "phi0" => phi0,
        "thetatgt" => thetatgt,
        "phitgt" => phitgt,
        "ef" => ef_norm,
        "rf" => rf_norm,
        "BAP" => guidance_logic[:type],
        "KBR" => KBR,
        "dsig" => 0.01 * π / 180,
        "sigf" => guidance_logic[:sigf],
        "KEF" => get(guidance_logic, :KEF, 0.0),
        "KLF" => get(guidance_logic, :KLF, 0.0),
        "vehicle" => vehicle,
    )

    if state.first_call
        delpsi0 = cal_delpsi(X0_norm, thetatgt, phitgt)
        state.BR = sign(delpsi0) > 0 ? -1 : 1
        state.first_call = false
    end

    if tc <= GAT / TU
        return 0.0
    end

    while true
        G = cal_z(X0_norm, state.sigs, ec_norm, auxdata)
        if abs(G) <= 1e3 / DU
            break
        end
        J = cal_J(X0_norm, state.sigs, ec_norm, auxdata)
        println("   Correction at t = ", round(tc * TU), " s")
        if abs(J) < 1e-6
            break
        end
        state.sigs -= G / J
    end

    siglmt = bank_constraints[:siglmt]

    if bank_constraints[:type] == 1 && abs(state.sigs) >= siglmt
        state.sigs = siglmt
    end

    BRprev = state.BR
    state.BR = cal_BR_prdt(X0_norm, ec_norm, state.sigs, BRprev, auxdata)


    if state.BR != BRprev
        Rgo = DU * 1e-3 * cal_sphdis(X0_norm[2], X0_norm[3], thetatgt, phitgt)
        println("       Bank reversal at Rgo = ", round(Rgo), " km")
    end

    bankcmd = state.sigs * state.BR

    if bank_constraints[:type] == 2
        sigdlmt = bank_constraints[:sigdlmt]
        sigddlmt = bank_constraints[:sigddlmt]
        tstp = 1.0
        minbnk = max(state.sigcmdprv - sigdlmt * tstp, 2 * state.sigcmdprv - state.sigcmdprv2 - sigddlmt * tstp^2)
        maxbnk = min(state.sigcmdprv + sigdlmt * tstp, 2 * state.sigcmdprv - state.sigcmdprv2 + sigddlmt * tstp^2)
        bankcmd = clamp(bankcmd, minbnk, maxbnk)

        if abs(bankcmd) >= siglmt
            bankcmd = siglmt
        end
    end

    state.sigcmdprv2 = state.sigcmdprv
    state.sigcmdprv = bankcmd

    return bankcmd
end


# Function for 3-degree of freedom equations of motion of the entry vehicle
# Input:
#   t       = time (independent variable)
#   X       = state vector [r, theta, phi, V, gamma, psi]
#   sigma   = bank angle
#   auxdata = dictionary containing various parameters for the simulation
# Output:
#   Xdot    = time derivative of the state vector

function aux_dyn!(dX, X, p, t)
    sigma, auxdata = p

    # Extract auxdata
    rp = auxdata["rp"]
    mu = auxdata["mu"]
    w = auxdata["w"]
    S = auxdata["S"]
    m = auxdata["m"]
    CL = auxdata["CL"]
    CD = auxdata["CD"]
    DU = auxdata["DU"]
    VU = auxdata["VU"]
    TU = auxdata["TU"]
    vehicle = auxdata["vehicle"]

    # State variables
    r = X[1]
    theta = X[2]
    phi = X[3]
    V = X[4]
    gamma = X[5]
    psi = X[6]

    # Aerodynamic forces and gravity
    Vreal = V * VU # m/s
    hreal = (r - rp) * DU # m
    # rho = cal_airdens(vehicle) # kg/m^3
    rho = cal_airdens(hreal, vehicle)  # kg/m^3
    q = 0.5 * rho[1] * Vreal^2
    D = q * S * CD / m
    L = q * S * CL / m
    D = D / DU * TU^2 # Unitless drag
    L = L / DU * TU^2 # Unitless lift
    g = mu / r^2     # Unitless gravity

    # Time domain 3-DoF entry dynamics
    rdot = V * sin(gamma)
    thetadot = V * cos(gamma) * sin(psi) / (r * cos(phi))
    phidot = (V / r) * cos(gamma) * cos(psi)
    Vdot = -D - g * sin(gamma) + w^2 * r * cos(phi) * (sin(gamma) * cos(phi) - cos(gamma) * sin(phi) * cos(psi))
    gamdot = (L / V) * cos(sigma) + (V / r - g / V) * cos(gamma) + 2 * w * cos(phi) * sin(psi) + w^2 * r / V * cos(phi) * (cos(gamma) * cos(phi) + sin(gamma) * sin(phi) * cos(psi))
    psidot = (L * sin(sigma)) / (V * cos(gamma)) + (V / r) * sin(psi) * cos(gamma) * tan(phi) - 2 * w * (tan(gamma) * cos(phi) * cos(psi) - sin(phi)) + w^2 * r / (V * cos(gamma)) * sin(phi) * cos(phi) * sin(psi)
    # Return the derivative of the state vector
    dX .= [rdot, thetadot, phidot, Vdot, gamdot, psidot]
    return dX
end

# Function for 2-degree of freedom equations of motion of the entry vehicle
# Input:
#   e       = energy-like variable
#   X       = state vector [s, r, gamma]
#   sigs    = initial bank angles
#   e0      = initial energy
#   auxdata = dictionary containing various parameters for the simulation
# Output:
#   Xdot    = time derivative of the state vector
function aux_dynE_2D!(dX, X, p, e)
    e0, sigs, auxdata = p

    # Extract auxdata
    rp = auxdata["rp"]
    mu = auxdata["mu"]
    S = auxdata["S"]
    m = auxdata["m"]
    CL = auxdata["CL"]
    CD = auxdata["CD"]
    DU = auxdata["DU"]
    VU = auxdata["VU"]
    TU = auxdata["TU"]
    ef = auxdata["ef"]
    BAP = auxdata["BAP"]
    vehicle = auxdata["vehicle"]

    # State variables
    s = X[1]
    r = X[2]
    gamma = X[3]
    V = sqrt(2 * (mu / r - e))

    # Aerodynamic forces and gravity
    Vreal = V * VU  # m/s
    hreal = (r - rp) * DU  # m
    rho = cal_airdens(hreal, vehicle)  # kg/m^3
    q = 0.5 * rho[1] * Vreal^2
    D = q * S * CD / m
    L = q * S * CL / m
    D = D / DU * TU^2  # Unitless drag
    L = L / DU * TU^2  # Unitless lift
    g = mu / r^2       # Unitless gravity

    # Bank angle parameterization
    if BAP == 1  # linear
        sigf = auxdata["sigf"]
        sig0 = sigs[1]
        sigma = sig0 + (sigf - sig0) / (ef - e0) * (e - e0)

    elseif BAP == 2  # exponential 1
        KEF = auxdata["KEF"]
        sigf = auxdata["sigf"]
        sig0 = sigs[1]
        sigma = exp(-(e - e0) / (ef - e0) * KEF) * (sig0 - sigf) + sigf

    elseif BAP == 3  # exponential 2
        KEF = auxdata["KEF"]
        sig0 = sigs[1]
        sigma = sig0 * exp(-(e - e0) / (ef - e0) * KEF)

    elseif BAP == 4  # logistic
        KLF = auxdata["KLF"]
        sig0 = sigs[1]
        sigma = 2 * sig0 / (1 + exp((e - e0) / (ef - e0) * KLF))
    end

    # Energy domain 2-DoF entry dynamics
    dedt = D * V
    dsdt = V * cos(gamma) / r
    drdt = V * sin(gamma)
    dgdt = (L / V) * cos(sigma) + (V / r - g / V) * cos(gamma)

    dsde = dsdt / dedt
    drde = drdt / dedt
    dgde = dgdt / dedt

    # Return the derivative of the state vector
    dX .= [dsde, drde, dgde]

    return dX
end

# Function for 3-degree of freedom equations of motion of the entry vehicle
# Input:
#   e       = energy-like variable
#   X       = state vector [r, theta, phi, V, gamma, psi]
#   sigs    = initial bank angles
#   BR      = bank angle reversal factor
#   e0      = initial energy
#   auxdata = dictionary containing various parameters for the simulation
# Output:
#   Xdot    = energy derivative of the state vector

function aux_dynE!(dX, X, p, e)
    e0, sigs, BR, auxdata = p

    # Extract auxdata
    rp = auxdata["rp"]
    mu = auxdata["mu"]
    w = auxdata["w"]
    S = auxdata["S"]
    m = auxdata["m"]
    CL = auxdata["CL"]
    CD = auxdata["CD"]
    DU = auxdata["DU"]
    VU = auxdata["VU"]
    TU = auxdata["TU"]
    ef = auxdata["ef"]
    BAP = auxdata["BAP"]
    vehicle = auxdata["vehicle"]

    # State variables
    r = X[1]
    theta = X[2]
    phi = X[3]
    V = X[4]
    gamma = X[5]
    psi = X[6]

    # Aerodynamic forces and gravity
    Vreal = V * VU  # m/s
    hreal = (r - rp) * DU  # m
    # rho = cal_airdens(vehicle) # kg/m^3
    rho = cal_airdens(hreal, vehicle)  # kg/m^3
    q = 0.5 * rho[1] * Vreal^2
    D = q * S * CD / m
    L = q * S * CL / m
    D = D / DU * TU^2  # Unitless drag
    L = L / DU * TU^2  # Unitless lift
    g = mu / r^2       # Unitless gravity

    # Bank angle parameterization
    if BAP == 1  # linear
        sigf = auxdata["sigf"]
        sig0 = sigs[1]
        sigma = sig0 + (sigf - sig0) / (ef - e0) * (e - e0)

    elseif BAP == 2  # exponential 1
        KEF = auxdata["KEF"]
        sigf = auxdata["sigf"]
        sig0 = sigs[1]
        sigma = exp(-(e - e0) / (ef - e0) * KEF) * (sig0 - sigf) + sigf

    elseif BAP == 3  # exponential 2
        KEF = auxdata["KEF"]
        sig0 = sigs[1]
        sigma = sig0 * exp(-(e - e0) / (ef - e0) * KEF)

    elseif BAP == 4  # logistic
        KLF = auxdata["KLF"]
        sig0 = sigs[1]
        sigma = 2 * sig0 / (1 + exp((e - e0) / (ef - e0) * KLF))
    end

    # Bank sign
    sigma *= BR

    # Energy domain 3-DoF entry dynamics
    dedt = D * V
    rdot = V * sin(gamma)
    thetadot = V * cos(gamma) * sin(psi) / (r * cos(phi))
    phidot = (V / r) * cos(gamma) * cos(psi)
    Vdot = -D - g * sin(gamma) + w^2 * r * cos(phi) * (sin(gamma) * cos(phi) - cos(gamma) * sin(phi) * cos(psi))
    gamdot = (L / V) * cos(sigma) + (V / r - g / V) * cos(gamma) + 2 * w * cos(phi) * sin(psi) + w^2 * r / V * cos(phi) * (cos(gamma) * cos(phi) + sin(gamma) * sin(phi) * cos(psi))
    psidot = (L * sin(sigma)) / (V * cos(gamma)) + (V / r) * sin(psi) * cos(gamma) * tan(phi) - 2 * w * (tan(gamma) * cos(phi) * cos(psi) - sin(phi)) + w^2 * r / (V * cos(gamma)) * sin(phi) * cos(phi) * sin(psi)

    # Return the energy derivative of the state vector
    dX .= [rdot, thetadot, phidot, Vdot, gamdot, psidot] ./ dedt

    return dX
end

# Function to calculate air density
# Input: 
#   h  = altitude in meters
#   pn = planet (1 for Earth, 2 for Mars)
# Output:
#   rho = air density in kg/m^3

# function cal_airdens(vehicle)
#     return Genesis.density(atmosphere(vehicle))
# end

function cal_airdens(h, vehicle_in)
    r_eq = equatorial_radius(planet(vehicle_in))
    currentPosition = planetodetic_altitude(vehicle_in) + r_eq
    currentLongitude = mod(longitude(vehicle_in) + 2pi, 2pi)
    currentLatitude = planetodetic_latitude(vehicle_in)

    set_position!(vehicle_in,
        radius=h + r_eq,
        longitude=currentLongitude,
        declination=currentLatitude)

    # Get the density at new position
    rho = Genesis.density(atmosphere(vehicle_in))

    # Revert to original position
    set_position!(vehicle_in,
        radius=currentPosition,
        longitude=currentLongitude,
        declination=currentLatitude)

    return rho
end

# ----------------------------------------------------------------------- #
# predictive lateral guidance based on
# Smith, K. M. (2016, February).
# Predictive lateral logic for numerical entry guidance algorithms.
# In AAS/AIAA Space Flight Mechanics Meeting (No. JSC-CN-35110-1).
# ----------------------------------------------------------------------- #
# input  = current state, previous bank sign,
#          parameterzied bank angle of the NPCG algorithm
# output = bank sign
function cal_BR_prdt(X0, e0, bank, BRprev, auxdata)

    # Extract KBR from auxdata
    KBR = auxdata["KBR"]

    # Crossrange error based on the current bank profile (mag + sign)
    _, CR1 = cal_z_3d(X0, e0, bank, BRprev, auxdata)

    # Crossrange error based on the opposite sign of bank profile
    _, CR2 = cal_z_3d(X0, e0, bank, -BRprev, auxdata)

    # Calculate the ratio of crossrange errors
    ratio = abs(CR1 / CR2)

    # Decide on the bank sign based on the ratio
    if ratio < KBR
        BR = BRprev
    else
        BR = -BRprev
    end

    return BR
end

# Bank reversal algorithm based on Tu et al. (2000)
# Input:
#   X        = current state vector
#   BR       = previous bank sign
#   case_num = case number for selecting guidance strategy
#   auxdata  = auxiliary data dict containing planetary parameters
# Output:
#   BR = updated bank sign

function cal_BR(X, BR, case_num, auxdata)

    # Extract auxiliary data
    rp = auxdata["rp"]
    DU = auxdata["DU"]
    thetatgt = auxdata["thetatgt"]
    phitgt = auxdata["phitgt"]
    dlpsT = auxdata["dlpsT"]

    # State variables
    r = X[1]
    theta = X[2]
    phi = X[3]
    V = X[4]
    gamma = X[5]
    psi = X[6]

    # sin/cos of angles
    cgam = cos(gamma)
    sgam = sin(gamma)
    cpsi = cos(psi)
    spsi = sin(psi)
    cth = cos(theta)
    sth = sin(theta)
    cphi = cos(phi)
    sphi = sin(phi)

    # ----------------------------------------------------------------------- #
    #          velocity vector projection onto local horizontal plane
    # ----------------------------------------------------------------------- #
    # velocity vector projection onto NE plane --> NED to ECEF
    Vp = V * cgam * [-cth * sphi * cpsi - sth * spsi;
             -sth * sphi * cpsi + cth * spsi;
             cphi * cpsi]
    UVp = Vp / norm(Vp)

    # position vector in inertial frame
    rnav = [r * cos(theta) * cphi;
        r * sin(theta) * cphi;
        r * sphi]

    # target vector in inertial frame
    rtgt = rp * [cos(thetatgt) * cos(phitgt);
        sin(thetatgt) * cos(phitgt);
        sin(phitgt)]

    Del = rtgt - rnav
    Ur = rnav / norm(rnav)

    # ----------------------------------------------------------------------- #
    #       target pointing vector projection onto local horizontal plane
    # ----------------------------------------------------------------------- #
    Delp = cross(cross(Ur, Del), Ur)
    UDelp = Delp / norm(Delp)

    # heading angle error magnitude
    delpsi = acos(dot(UVp, UDelp))

    # Determine whether target is on the left or right
    sign_ck = dot(cross(UVp, UDelp), Ur)

    # current range-to-go
    R2go = cal_sphdis(theta, phi, thetatgt, phitgt) * DU * 1e-3  # km

    # Define heading angle error thresholds based on case_num
    if case_num == 1  # Mars robotic
        del_psi_lim = R2go > 100 ? dlpsT[1] : dlpsT[2]
    elseif case_num == 2 || case_num == 3  # Mars manned / Mars manned new
        del_psi_lim = R2go > 150 ? dlpsT[1] : dlpsT[2]
    elseif case_num == 4  # Apollo 10
        del_psi_lim = R2go > 150 ? dlpsT[1] : dlpsT[2]
    elseif case_num == 5  # CAV-H
        del_psi_lim = R2go > 300 ? dlpsT[1] : dlpsT[2]
    end

    # If current heading angle error exceeds the limit, adjust bank sign
    if abs(delpsi) >= del_psi_lim
        if sign_ck > 0
            # Target is on the left side of the flight direction -> turn left -> negative roll
            BR = -1
        else
            # Target is on the right side of the flight direction -> turn right -> positive roll
            delpsi = -abs(delpsi)
            BR = 1
        end
    end

    return BR
end

# Function to calculate heading angle error (del_psi)
# Input:
#   X      = state vector of vehicle [r, theta, phi, V, gamma, psi]
#   thetat = target latitude in radians
#   phit   = target longitude in radians
# Output:
#   del_psi = heading angle error in radians

function cal_delpsi(X::Vector{Float64}, thetat::Float64, phit::Float64)

    # State variables
    r = X[1]
    theta = X[2]
    phi = X[3]
    V = X[4]
    gamma = X[5]
    psi = X[6]

    # sin/cos of angles
    cgam = cos(gamma)
    sgam = sin(gamma)
    cpsi = cos(psi)
    spsi = sin(psi)
    cth = cos(theta)
    sth = sin(theta)
    cphi = cos(phi)
    sphi = sin(phi)
    ctht = cos(thetat)
    stht = sin(thetat)
    cphit = cos(phit)
    sphit = sin(phit)

    # ------ projection of velocity vector onto local horizontal plane ------ #
    # 1. Get velocity vector in MCR component (NED frame)
    # 2. Projection onto the local-horizontal frame (= North East plane)
    #    Third component = 0
    # 3. Express the velocity vector in planet-fixed frame (ex. ECEF)
    # 4. Normalize the velocity vector
    Vp = V * cgam * [-cth * sphi * cpsi - sth * spsi;
             -sth * sphi * cpsi + cth * spsi;
             cphi * cpsi]
    UVp = Vp / norm(Vp)

    # Position vector
    r_nav = [cth * cphi; sth * cphi; sphi]
    Ur = r_nav / norm(r_nav)

    # Target vector on the ground
    r_tgt = [ctht * cphit; stht * cphit; sphit]
    Del = r_tgt - r_nav

    # Project Del onto the horizontal plane
    Delp = cross(cross(Ur, Del), Ur)
    UDelp = Delp / norm(Delp)

    # Heading error magnitude
    del_psi = acos(dot(UVp, UDelp))

    # Heading error direction
    # If Vp is on the right side of Delp
    # --> cross(UVp, UDelp) // Ur   --> + dot product
    # If Vp is on the left side of Delp
    # --> cross(UVp, UDelp) // -Ur   --> - dot product
    sign_ck = sign(dot(cross(UVp, UDelp), Ur))

    if sign_ck < 0
        del_psi = -abs(del_psi)
    end

    return real(del_psi)
end

# Function to calculate downrange (DR) and crossrange (CR)
# Input:
#   theta0, phi0 = longitude and latitude of initial location (in radians)
#   thetaf, phif = longitude and latitude of target location (in radians)
#   theta, phi   = longitude and latitude of current location (in radians)
# Output:
#   DR = downrange in radians
#   CR = crossrange in radians

function cal_drdc(theta0::Float64, phi0::Float64, thetaf::Float64, phif::Float64, theta::Float64, phi::Float64)

    # sin/cos of angles
    cth0 = cos(theta0)
    sth0 = sin(theta0)
    cph0 = cos(phi0)
    sph0 = sin(phi0)
    cthf = cos(thetaf)
    sthf = sin(thetaf)
    cphf = cos(phif)
    sphf = sin(phif)
    cth = cos(theta)
    sth = sin(theta)
    cph = cos(phi)
    sph = sin(phi)

    # Calculate unit distance from X0 to X on the unit sphere
    c = cal_sphdis(theta0, phi0, theta, phi)

    # Position vectors
    X0 = [cth0 * cph0; sth0 * cph0; sph0]  # Initial position vector
    Xt = [cthf * cphf; sthf * cphf; sphf]  # Target position vector
    X = [cth * cph; sth * cph; sph]       # Current position vector

    # Normal vector of X0 and Xt
    nX0Xt = cross(X0, Xt)
    nX0Xt /= norm(nX0Xt)

    # Normal vector of X0 and X
    if theta0 == theta && phi0 == phi  # If at initial point
        nX0X = nX0Xt
    else
        nX0X = cross(X0, X)
        nX0X /= norm(nX0X)
    end

    # Included angle between X0toX and X0toXt
    A = acos(dot(nX0X, nX0Xt))
    A = real(A)

    # Crossrange (CR): sin(a) = sin(A) * sin(c)
    a = asin(sin(A) * sin(c))
    CR = a

    # Downrange (DR): cos(c) = cos(a) * cos(b)
    b = acos(cos(c) / cos(a))
    DR = b

    # Define crossrange sign
    temp0 = cross(nX0X, nX0Xt)
    temp1 = temp0 / norm(temp0)
    temp2 = dot(temp1, X0)

    if sign(temp2) < 0
        CR = -CR
    end

    return DR, CR
end

# ----------------------------------------------------------------------- #
# calculate derivative using the parameterized bank angle
# ----------------------------------------------------------------------- #
# input  = current state, guidance
# output = Jacobian

function cal_J(X0, sigs, e0, auxdata)

    # Extract dsig from auxdata
    dsig = auxdata["dsig"]

    # Small deviation in control variable
    u1 = sigs[1] + dsig
    G1 = cal_z(X0, u1, e0, auxdata)

    u2 = sigs[1] - dsig
    G2 = cal_z(X0, u2, e0, auxdata)

    # Calculate the Jacobian
    J = (G1[1] - G2[1]) / (2 * dsig)

    return J
end

# Function to transform coordinates from PCI (planet-centered inertial) 
# to PCR (planet-centered rotating) based on the NED frame.
# Input:
#   ri, thetai, phii = radial distance, latitude, and longitude in PCI
#   Vi, gammai, psii = velocity, flight path angle, and heading angle in PCI
#   w                = planet rotation rate
# Output:
#   rn, thetan, phin = radial distance, latitude, and longitude in PCR
#   Vn, gamman, psin = velocity, flight path angle, and heading angle in PCR

function cal_pci2pcr(ri::Float64, thetai::Float64, phii::Float64, Vi::Float64, gammai::Float64, psii::Float64, w::Float64)

    # Angles
    sth = sin(thetai)
    cth = cos(thetai)
    sphi = sin(phii)
    cphi = cos(phii)
    cgam = cos(gammai)
    sgam = sin(gammai)
    cpsi = cos(psii)
    spsi = sin(psii)

    # Assumption: PCI and PCR are overlapped at the moment
    thetan = thetai
    phin = phii

    # Radial vector in PCI
    rI = [ri * cphi * cth;
        ri * cphi * sth;
        ri * sphi]

    # Velocity vector in PCI
    vI = [Vi * (cth * cphi * sgam - sth * spsi * cgam - cth * sphi * cgam * cpsi);
        Vi * (sth * cphi * sgam + cth * spsi * cgam - sth * sphi * cgam * cpsi);
        Vi * (sphi * sgam + cphi * cgam * cpsi)]

    # Planet rotation vector in PCI
    wI = [0.0; 0.0; w]

    # Velocity vector in PCR (relative velocity)
    vR = vI - cross(wI, rI)

    # Radius
    rn = norm(rI)

    # Relative speed
    Vn = norm(vR)

    # Planet-relative flight-path angle
    gamman = asin(dot(rI, vR) / (rn * Vn))

    # Planet-relative heading angle
    sgamn = sin(gamman)
    cgamn = cos(gamman)
    psin = acos((vR[3] - Vn * sgamn * sphi) / (Vn * cgamn * cphi))

    return rn, thetan, phin, Vn, gamman, psin
end

# Function to calculate great-circle distance on a unit spherical surface
# Input:
#   theta1, phi1 = longitude and latitude of the first location (in radians)
#   theta2, phi2 = longitude and latitude of the second location (in radians)
# Output:
#   ang = angle distance in radians

function cal_sphdis(theta1::Float64, phi1::Float64, theta2::Float64, phi2::Float64)
    ang = acos(sin(phi2) * sin(phi1) + cos(phi2) * cos(phi1) * cos(theta2 - theta1))
    return ang
end

# ----------------------------------------------------------------------- #
# calculate downrange and crossrange at the final energy state
# using the parameterized bank angle
# ----------------------------------------------------------------------- #
# input  = current state
# output = downrange and crossrange

# calculate downrange (DR) and crossrange (CR) at the final energy state
function cal_z_3d(X0, e0, convar, BR, auxdata)

    # Extract auxdata from the dictionary
    theta0 = auxdata["theta0"]
    phi0 = auxdata["phi0"]
    thetatgt = auxdata["thetatgt"]
    phitgt = auxdata["phitgt"]
    ef = auxdata["ef"]

    # Energy span
    espan = [e0, ef]
    # parameters
    p = [e0, convar, BR, auxdata]

    prob = ODEProblem(aux_dynE!, X0, espan, p)

    # Solve the ODE
    sol = solve(prob, BS5(), reltol=1e-6, abstol=1e-9)

    # Extract time points and solution
    X = sol.u  # Solution at those time points    

    # Extract the final theta and phi values from the solution
    theta = X[end][2]
    phi = X[end][3]

    # Calculate downrange (DR) and crossrange (CR)
    DR, CR = cal_drdc(theta0, phi0, thetatgt, phitgt, theta, phi)

    return DR, CR
end

# calculate range-to-go using the parameterized bank angle
function cal_z(X0, sigs, e0, auxdata)
    # Extract auxdata from the dictionary
    thetatgt = auxdata["thetatgt"]
    phitgt = auxdata["phitgt"]
    ef = auxdata["ef"]

    # State vector
    r0 = X0[1]
    theta0 = X0[2]
    phi0 = X0[3]
    gamma0 = X0[5]

    # Current range-to-go
    sf = cal_sphdis(theta0, phi0, thetatgt, phitgt)

    # Initial condition
    x0 = [0.0, r0, gamma0]

    # Energy span
    ans = abs(e0 - ef) / 10000
    espan = [e0, ef]

    # parameters
    p = [e0, sigs, auxdata]

    prob = ODEProblem(aux_dynE_2D!, x0, espan, p)

    # Solve the ODE
    sol = solve(prob, BS5(), saveat=ans, reltol=1e-6, abstol=1e-9)
    # Extract time points and solution
    X = sol.u  # Solution at those time points    

    # sef = X(end,1) - sf
    z = X[end][1] - sf

    return z
end

end # module
