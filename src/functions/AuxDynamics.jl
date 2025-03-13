module AuxDynamics

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
    pn = auxdata["pn"]
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
    rho = cal_airdens(hreal, pn) # kg/m^3
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
    pn = auxdata["pn"]
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

    # State variables
    s = X[1]
    r = X[2]
    gamma = X[3]
    V = sqrt(2 * (mu / r - e))

    # Aerodynamic forces and gravity
    Vreal = V * VU  # m/s
    hreal = (r - rp) * DU  # m
    rho = cal_airdens(hreal, pn)  # kg/m^3
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
    pn = auxdata["pn"]
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
    rho = cal_airdens(hreal, pn)  # kg/m^3
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

end