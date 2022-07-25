using Formatting
include("detect_ranges.jl")


ps = [ 0, 0.5, 1, 1.5 , 2, 2.5]
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

for (i, p) in enumerate(ps)

    coe = (
        c = 3.6e3,
        μ = 3.0,
        ν = 1.0,
        p = p,
    )

    ψs = collect(range(ψ_rng..., length=1001))
    ξs = coe.μ / coe.ν * ( (1 .- ψs / coe.μ) / coe.p .- 1 ./ (1 .+ abs.(ψs)) )

    d_dydτ_dy = - ( 1 .+ abs.(ψs) ) + coe.μ * sign.(ψs) .* (1 .- (ψs .- coe.ν * ps * ξ) / coe.μ )
    stable = d_dydτ_dy .< 0

    vals, rngs = detectRanges(stable)
    
    color = ["black", "red", "blue", "green", "orange", "purple"][i]
    for (k, rng) in enumerate(rngs)

        args = Dict()
        if k==1
            args[:label] = format("\$p = {:f} \$", coe.p)
        else
            args[:label] = nothing
        end

        args[:linestyle] = (vals[k] == 0) ? "dotted" : "solid"
        args[:color] = color
        ax.plot(ξs[rng], ψs[rng]; args...)

    end
        
    #ax.text(label_positions[i]..., format("\$\\xi = {:.1f} \$", ξ), va="center", ha="center", color=color, size=15)

end

ax.set_xlim([-10, 10])
ax.set_ylim([-1, 5])
ax.legend(loc="center right",)
fig.savefig("figures/figure-stommel_bifurcation_analytical_xi.png", dpi=300)

plt.show()
