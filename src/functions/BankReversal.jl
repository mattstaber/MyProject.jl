module BankReversal

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

end