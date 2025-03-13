module calAirdens

# Function to calculate air density
# Input: 
#   h  = altitude in meters
#   pn = planet (1 for Earth, 2 for Mars)
# Output:
#   rho = air density in kg/m^3

function cal_airdens(h, pn)
    rho = zeros(Float64, length(h))  # Initialize a vector for densities

    if pn == 1
        # Earth atmosphere (Hicks 2009)
        beta = 0.14  # 1/km
        rhos = 1.225
        rho .= rhos .* exp.(-beta .* h .* 1e-3)  # Use broadcasting (.)

    elseif pn == 2
        # Mars atmosphere (AAS 18-485)
        beta1 = 559.351005946503
        beta2 = 188.95110711075

        T = (1.4e-13) .* (h .^ 3) .- (8.85e-9) .* (h .^ 2) .- (1.245e-3) .* h .+ 205.3645
        rho0 = beta1 ./ (beta2 .* T)
        beta = -0.000105 .* h
        rho .= rho0 .* exp.(beta)  # Use broadcasting (.)

    else
        error("Invalid planet number. Use 1 for Earth or 2 for Mars.")
    end

    return rho
end

end