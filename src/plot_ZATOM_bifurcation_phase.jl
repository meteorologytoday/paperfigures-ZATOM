using Formatting
using Roots

# ξ  γ_left γ_right
typeA_pts = [
    -0.50  0.2321  0.4585 ;
    -0.25  0.1817  0.4048 ;
     0.00  0.1451  0.3396 ;
     0.25  0.1125  0.2539 ;
     0.50  0.0893  0.1508 ;
     0.75  0.0726  0.0924 ;
     1.00  0.0608  0.0672 ;
]

typeB_pts = [
     0.95  0.09008  0.14209 ;
     1.00  0.08407  0.11485 ;
     1.05  0.07893  0.09770 ;
     1.10  0.07598  0.08627 ;
     1.15  0.07112  0.07750 ;
]


println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel("\$\\gamma\$ [ Sv ]", fontsize=25)
ax.set_ylabel("\$\\xi\$", fontsize=25)
ax.grid()
#ax.set_ylim([-1, 6.2])
#ax.set_xlim([0.0, 6])
ax.invert_yaxis()


ax.fill_betweenx(typeA_pts[:, 1], typeA_pts[:, 2], typeA_pts[:, 3], facecolor="lightsalmon", edgecolor="orangered", alpha=0.8, linewidth=1, zorder=10, label="\$\\mu_A\$")
ax.fill_betweenx(typeB_pts[:, 1], typeB_pts[:, 2], typeB_pts[:, 3], facecolor="dodgerblue", edgecolor="blue", alpha=0.8, linewidth=1, zorder=10, label="\$\\mu_B\$")

ax.legend(loc="lower right")
ax.set_title("(b) Multiple equilibria phase diagram of ZATOM")
fig.savefig("figures/figure-ZATOM_bifurcation_phase.png", dpi=300)

plt.show()
