using Formatting
using Roots

# ξ  γ_left γ_right
typeA_pts = [
     0.95  0.0901   0.1422;
     1.00  0.0840   0.1147;
     1.05  0.0789   0.0980;
     1.10  0.0758   0.0865;
     1.15  0.0711   0.0774;
     1.20  0.0676   0.0711;
     1.25  0.0638   0.0662;
     1.30  0.0600   0.0619;
     1.35  0.0570   0.0583;
     1.40  0.0542   0.0549;
     1.45  0.0516   0.0522;
     1.50  0.0488   0.0495;
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
#ax.fill_betweenx(typeB_pts[:, 1], typeB_pts[:, 2], typeB_pts[:, 3], facecolor="dodgerblue", edgecolor="blue", alpha=0.8, linewidth=1, zorder=10, label="\$\\mu_B\$")

#ax.legend(loc="lower right")
ax.set_title("Multiple equilibria phase diagram of ZATOM")
fig.savefig("figures/figure-ZATOM_bifurcation_phase.png", dpi=300)

plt.show()
