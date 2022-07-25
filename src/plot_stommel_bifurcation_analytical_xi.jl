using Formatting
include("detect_ranges.jl")


ps = [ 0.1, 0.5, 1.0, 1.5 , 2, 2.5, 1e5]
ψ0 = 0.0
ψ_rng = [-4, 6]

μ = 3.0
ν = 1.0

label_positions = [
    (-5.00, 2.75),
    (-4.00, 0.4),
    (-1.30, 0.4),
    nothing,
    nothing,
    ( 2.20, 0.4),
    ( 2.00, 1.60),
]


println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)



ax.set_xlabel("\$\\xi\$", fontsize=25)
ax.set_ylabel("\$\\Psi\$", fontsize=25)
ax.grid()

for (i, p) in enumerate(ps)

    coe = (
        c = 3.6e3,
        μ = μ,
        ν = ν,
        p = p,
    )

    ψs = collect(range(ψ_rng..., length=1001))
    ξs = coe.μ / coe.ν * ( 1 ./ (1 .+ abs.(ψs)) .- (1 .- ψs / coe.μ) / coe.p )

    d_dydτ_dy = - ( 1 .+ abs.(ψs) ) + coe.μ * sign.(ψs) .* (1 .- (ψs .- coe.ν .* coe.p .* ξs) / coe.μ )
    stable = d_dydτ_dy .< 0

    vals, rngs = detectRanges(stable)
    
    color = ["black", "red", "orange", "green", "blue", "purple", "gray"][i]
    for (k, rng) in enumerate(rngs)

        args = Dict()
        if k==1
            #args[:label] = format("\$p = {:f} \$", coe.p)
        else
            args[:label] = nothing
        end

        args[:linestyle] = (vals[k] == 0) ? "dotted" : "solid"
        args[:color] = color
        ax.plot(ξs[rng], ψs[rng]; args...)

    end
    
    if label_positions[i] != nothing

        if p < 1e3 
            txt = format("\$p = {:.1f} \$", p)
        else
            txt = format("\$p \\rightarrow \\infty \$", p)
        end
        ax.text(label_positions[i]..., txt, va="center", ha="center", color=color, size=15)
    end

end

# ξ, ψ
C1 = ( μ / (1 + μ) / ν , μ )
ax.scatter(C1..., s=20, marker="o", color="black", zorder=99)
ax.annotate("\$ C_1 = \\left( \\frac{\\mu}{\\left( 1 + \\mu \\right) \\nu }, \\mu \\right) \$", xy=C1,  xycoords="data",
            xytext=(-2, 3.5), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=13,
)

_p = 0.5
C2 = ( μ / ν * (1 - 1 / _p), 0 )
ax.scatter(C2..., s=20, marker="o", color="black", zorder=99)
ax.annotate("\$ C_2 = \\left( \\frac{\\mu}{\\nu} \\left( 1 - \\frac{1}{p} \\right), 0 \\right) \$", xy=C2,  xycoords="data",
            xytext=(-3, -.7), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=13,
)
ax.set_xlim([-6, 4])
ax.set_ylim([-1, 5])
ax.text(0.05, 0.95, "(b)", size=25, va="top", ha="left", transform=ax.transAxes)
#ax.legend(loc="center right",)

fig.savefig("figures/figure-stommel_bifurcation_analytical_xi.png", dpi=300)

plt.show()
