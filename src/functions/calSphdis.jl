module calSphdis

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

end