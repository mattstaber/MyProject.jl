module calPci2pcr

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

end