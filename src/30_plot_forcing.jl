using Formatting

include("ZATOM_calculation_setup.jl")



println("Loading PyPlot")
println("Setting backend as Agg...")
ENV["MPLBACKEND"]="agg"
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, figsize=(6,4), constrained_layout=true)

ax_twinx = ax.twinx()

ϕ = collect(range(ϕs, ϕn, length=N))
ϕ_hlat = ϕ[ϕ .>= ϕc]



ln_T_sfc = ax.plot(rad2deg.(ϕ), T_sfc.(ϕ), "k-")
ln_σ = ax_twinx.plot(rad2deg.(ϕ), PmE_ξ0.(ϕ) * factor_PmE_mstommyear,   "r-")
ln_σ_w = ax_twinx.plot(rad2deg.(ϕ_hlat), PmE_w.(ϕ_hlat) * factor_PmE_mstommyear, ls="-.", c="red")
ln_σ_e = ax_twinx.plot(rad2deg.(ϕ_hlat), PmE_e.(ϕ_hlat) * factor_PmE_mstommyear, ls="--", c="red")

ax.set_xlabel("\$\\phi\$ [ \${}^{\\circ}\$N ]")

ax.set_ylabel("\$ T_{\\mathrm{sfc}} \$ [ \${}^{\\circ}\\mathrm{C} \$ ]")
ax_twinx.set_ylabel("PmE [ \$ \\mathrm{m} / \\mathrm{s}  \$ ]", color="red")
ax_twinx.spines["right"].set_color("red")

ax_twinx.text(55, 75, "\$ \\mathrm{PmE}_e \$", size=12, color="red", ha="center", va="center")
ax_twinx.text(65, 25, "\$ \\mathrm{PmE}_w \$", size=12, color="red", ha="center", va="center")

ax.set_yticks([0, 5, 10, 15, 20, 25, 30, 35])
ax_twinx.set_yticks([-75, -50, -25, 0, 25, 50, 75, 100])
ax.set_ylim([0, 35])
ax_twinx.set_ylim([-75, 100])
ax_twinx.yaxis.label.set_color("red")
ax_twinx.tick_params(axis="y", colors="red")
ax.grid()

ax.set_title("(a)", size=20)

output_file = "figures/quantitative-forcing.svg"
println("Output file: $output_file")
fig.savefig(output_file, dpi=300)

plt.show()
