module calDrdc

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

end