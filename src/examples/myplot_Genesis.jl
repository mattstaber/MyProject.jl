function plot_edl_data(sigmadata, V, h, DR, theta, phi, time, gamma, psi, CR, A, q, qdot, r2d)
    # Create a figure
    f = Figure(size=(1200, 675))

    # Plot sigma vs time
    ax1 = Axis(f[1, 1], xlabel="t (s)", ylabel="sigma (deg)")
    lines!(ax1, sigmadata[:, 1], sigmadata[:, 2] * r2d, color=:blue, linewidth=2)

    # Plot V vs h
    ax2 = Axis(f[1, 2], xlabel="V (km/s)", ylabel="h (km)")
    lines!(ax2, V, h, color=:blue, linewidth=2)

    # Plot DR vs h
    ax3 = Axis(f[1, 3], xlabel="Downrange (km)", ylabel="h (km)")
    lines!(ax3, DR, h, color=:blue, linewidth=2)

    # Plot theta vs phi
    ax4 = Axis(f[1, 4], xlabel="theta (deg)", ylabel="phi (deg)")
    lines!(ax4, theta, phi, color=:blue, linewidth=2)

    # Plot h vs time
    ax5 = Axis(f[2, 1], xlabel="t (s)", ylabel="h (km)")
    lines!(ax5, time, h, color=:blue, linewidth=2)

    # Plot V vs time
    ax6 = Axis(f[2, 2], xlabel="t (s)", ylabel="V (km/s)")
    lines!(ax6, time, V, color=:blue, linewidth=2)

    # Plot gamma vs time
    ax7 = Axis(f[2, 3], xlabel="t (s)", ylabel="gamma (deg)")
    lines!(ax7, time, gamma, color=:blue, linewidth=2)

    # Plot psi vs time
    ax8 = Axis(f[2, 4], xlabel="t (s)", ylabel="psi (deg)")
    lines!(ax8, time, psi, color=:blue, linewidth=2)

    # Plot Crossrange vs time
    ax9 = Axis(f[3, 1], xlabel="t (s)", ylabel="Crossrange (km)")
    lines!(ax9, time, CR, color=:blue, linewidth=2)

    # Plot A vs time with limits
    ax10 = Axis(f[3, 2], xlabel="t (s)", ylabel="A (g)")
    lines!(ax10, time, A, color=:blue, linewidth=2)

    # Plot q vs time with limits
    ax11 = Axis(f[3, 3], xlabel="t (s)", ylabel="q (kPa)")
    lines!(ax11, time, q * 1e-3, color=:blue, linewidth=2)

    # Plot qdot vs time with limits
    ax12 = Axis(f[3, 4], xlabel="t (s)", ylabel="Qdot (kW/m^2)")
    lines!(ax12, time, qdot, color=:blue, linewidth=2)

    # Return the figure
    return f
end