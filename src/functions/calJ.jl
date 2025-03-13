module calJ

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

end