using Formatting
using Roots

# ξ  γ_left γ_right
fixed_ξ = [
     -0.95  0.0901   0.1422;
     -1.00  0.0840   0.1147;
     -1.05  0.0789   0.0980;
     -1.10  0.0758   0.0865;
     -1.15  0.0711   0.0774;
     -1.20  0.0676   0.0711;
     -1.25  0.0638   0.0662;
     -1.30  0.0600   0.0619;
     -1.35  0.0570   0.0583;
     -1.40  0.0542   0.0549;
     -1.45  0.0516   0.0522;
     -1.50  0.0488   0.0495;
]

fixed_γ = [
    0.04    -1.741 -1.737;
    0.05    -1.486 -1.478;
    0.06    -1.322 -1.296;
    0.07    -1.207 -1.167;
    0.08    -1.133 -1.050;
    0.09    -1.082 -0.949;
    0.10    -1.042 -0.876;
    0.11    -1.011 -0.821;
    0.12    -0.987 -0.775;
    0.13    -0.969 -0.724;
    0.14    -0.952 -0.684;
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
#ax.invert_yaxis()


#ax.fill_betweenx(fixed_ξ[:, 1], fixed_ξ[:, 2], fixed_ξ[:, 3], facecolor="lightsalmon", edgecolor="orangered", alpha=0.8, linewidth=1, zorder=10, label="Folding in \$\\xi\$")
#ax.fill_between(fixed_γ[:, 1], fixed_γ[:, 2], fixed_γ[:, 3], facecolor="dodgerblue",  edgecolor="blue",      alpha=0.8, linewidth=1, zorder=10, label="Folding in \$\\gamma\$")

ax.fill_betweenx(fixed_ξ[:, 1], fixed_ξ[:, 2], fixed_ξ[:, 3], hatch="..", facecolor="none",  edgecolor="blue", alpha=0.8, linewidth=1, zorder=10, label="Folding in \$\\xi\$")
ax.fill_between(fixed_γ[:, 1], fixed_γ[:, 2], fixed_γ[:, 3], hatch="//",  facecolor="none",  edgecolor="orangered",      alpha=0.8, linewidth=1, zorder=10, label="Folding in \$\\gamma\$")

#ax.fill_betweenx(typeB_pts[:, 1], typeB_pts[:, 2], typeB_pts[:, 3], facecolor="dodgerblue", edgecolor="blue", alpha=0.8, linewidth=1, zorder=10, label="\$\\mu_B\$")

ax.legend(loc="lower right")
ax.set_title("Multiple equilibria phase diagram of ZATOM")
fig.savefig("figures/figure-ZATOM_bifurcation_phase.png", dpi=300)

plt.show()
