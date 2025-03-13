module EntryAnalysis

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

end