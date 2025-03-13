module calDelpsi

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

end