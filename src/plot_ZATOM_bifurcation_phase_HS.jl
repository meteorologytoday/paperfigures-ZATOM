using Formatting
using Roots


regime = Dict(
    "standard" => Dict(
        
        "label" => "\$H_S = 60 \\mathrm{m}\$",
        "label_pos" => (0.1, -0.6),
        "colors" => Dict(
            "fixed_ξ" => "blue", 
            "fixed_γ" => "orangered", 
        ),

        # ξ  γ_left γ_right
        "fixed_ξ" => [
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
        ],

        "fixed_γ" => [
        #    0.04    -1.741 -1.737;
        #    0.05    -1.486 -1.478;
        #    0.06    -1.322 -1.296;
        #    0.07    -1.207 -1.167;
            0.08    -1.133 -1.050;
            0.09    -1.082 -0.949;
            0.10    -1.042 -0.876;
            0.11    -1.011 -0.821;
            0.12    -0.987 -0.775;
            0.13    -0.969 -0.724;
            0.14    -0.952 -0.684;
        ],
    ),

    "H295m" => Dict(

        "label" => "\$H_S = 295 \\mathrm{m}\$",
        "label_pos" => (0.5, -0.1),
        "colors" => Dict(
            "fixed_ξ" => "blue", 
            "fixed_γ" => "orangered", 
        ),

        # ξ  γ_left γ_right
        "fixed_ξ" => [
             0.5000      0.8323      0.8427      ;
             0.4000      0.7931      0.8049      ;
             0.3000      0.7568      0.7723      ;
             0.2000      0.7220      0.7531      ;
             0.1000      0.6909      0.7331      ;
             0.0000      0.6583      0.7124      ;
            -0.1000      0.6243      0.6894      ;
            -0.2000      0.6006      0.6650      ;
            -0.6000         NaN      0.5060      ;
            -0.7000      0.2220      0.3850      ;
            -0.8000      0.1910      0.2800      ;
            -0.9000      0.1620      0.2060      ;
            -1.0000      0.1400      0.1610      ;
            -1.1000      0.1210      0.1290      ;
            -1.2000      0.1060      0.1090      ;
            -1.3000      0.0930      0.0930      ;
        ],

        "fixed_γ" => [
            0.075      -1.487      -1.483      ;
            0.100      -1.250      -1.238      ;
            0.150      -1.028      -0.950      ;
            0.200      -0.911      -0.764      ;
            0.250      -0.837      -0.624      ;
            0.300      -0.776      -0.506      ;
            0.350      -0.730      -0.407      ;
            0.400      -0.687      -0.328      ;
            0.450      -0.643      -0.252      ;
            0.500      -0.605      -0.217      ;
#            0.550      -0.544      NaN         ;
            0.600      -0.429      -0.202      ;
            0.650      -0.258      -0.025      ;
            0.700      -0.053       0.131      ;
            0.750       0.188       0.285      ;
            0.775       0.316       0.353      ;
            0.825       0.456       0.486      ;
            0.850       0.525       0.542      ;
            0.950       0.758       0.770      ;
            0.975       0.809       0.817      ;
        ],
    ),
)




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

for (key, data) in regime

    println("Doing key = $key")

    fixed_ξ = data["fixed_ξ"]
    fixed_γ = data["fixed_γ"]
    label = data["label"]

    colors = data["colors"]
 
    ax.fill_betweenx(fixed_ξ[:, 1], fixed_ξ[:, 2], fixed_ξ[:, 3], hatch="..", facecolor="none",  edgecolor=colors["fixed_ξ"], alpha=0.8, linewidth=1, zorder=10, label="[$label] Folding along fixed \$\\xi\$")
    ax.fill_between(fixed_γ[:, 1], fixed_γ[:, 2], fixed_γ[:, 3], hatch="//",  facecolor="none",  edgecolor=colors["fixed_γ"], alpha=0.8, linewidth=1, zorder=10, label="[$label] Folding along fixed \$\\gamma\$")

    ax.text(data["label_pos"]..., data["label"], size=12, ha="center", va="center")
end

#ax.legend(loc="lower right")
fig.savefig("figures/figure-ZATOM_bifurcation_phase_HS.png", dpi=300)

plt.show()
