using Formatting
include("detect_ranges.jl")


ξ_θs = [
    ( 1.0, 0.0),
    ( 0.5, 0.0),
    ( 0.0, 0.0),
    (-5.0, 0.0),
]
ψ0 = 0.0
ψ_rng = [-4, 6]


label_positions = [
    (1.25, 3.7),
    (1.85, 2.6),
    (1.55, 1.7),
    (0.70, 0.5),
]


println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)



ax.set_xlabel("\$p\$", fontsize=25)
ax.set_ylabel("\$\\Psi\$", fontsize=25)
ax.grid()

for (i, (ξ, θ)) in enumerate(ξ_θs)

    coe = (
        c = 3.6e3,
        μ = 3.0,
        ν = 1.0,
        θ = θ,
        ξ = ξ,
    )

    ψs = collect(range(ψ_rng..., length=1001))
    ps = (1 .+ abs.(ψs)) .* (1 .- ψs / coe.μ) ./ ( 1 .- coe.ν * ξ / coe.μ * ( 1 .+ abs.(ψs)) )

    d_dydτ_dy = - ( 1 .+ abs.(ψs) ) + coe.μ .* (1 .+ coe.θ) * sign.(ψs) .* (1 .- (ψs .- coe.ν * ps * ξ) / (coe.μ .* (1 .+ coe.θ)))
    stable = d_dydτ_dy .< 0

    vals, rngs = detectRanges(stable)
    
    color = ["black", "red", "orange", "green", "blue", "purple"][i]
    for (k, rng) in enumerate(rngs)

        args = Dict()
        if k==1
            #args[:label] = format("\$\\left(\\xi, \\theta \\right) = \\left( {:.1f}, {} \\right) \$", ξ, θ)
            #args[:label] = format("\$\\xi = {:.1f} \$", ξ)
        else
            args[:label] = nothing
        end

        args[:linestyle] = (vals[k] == 0) ? "dotted" : "solid"
        args[:color] = color
        ax.plot(ps[rng], ψs[rng]; args...)

    end
        
    ax.text(label_positions[i]..., format("\$\\xi = {:.1f} \$", ξ), va="center", ha="center", color=color, size=15)

end

ax.set_xlim([0, 3])
ax.set_ylim([-1, 5])
ax.text(0.05, 0.95, "(a)", size=25, va="top", ha="left", transform=ax.transAxes)
#ax.legend(loc="center right",)
fig.savefig("figures/figure-stommel_bifurcation_analytical_p.png", dpi=300)

plt.show()
