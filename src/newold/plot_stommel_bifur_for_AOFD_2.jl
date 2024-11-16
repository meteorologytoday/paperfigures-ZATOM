using Formatting
include("detect_ranges.jl")


ξ_θs = [
    ( -1.0, 0.0),
    ( -0.5, 0.0),
    (  0.0, 0.0),
    (  5.0, 0.0),
]
ψ0 = 0.0
ψ_rng = [-4, 6]

println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)


ax.set_xlabel("\$\\gamma\$", fontsize=25)
ax.set_ylabel("\$\\psi\$", fontsize=25)
ax.grid()

cnt = 1
for (ξ, θ) in ξ_θs

    coe = (
        c = 3.6e3,
        μ = 3.0,
        ν = 1.0,
        θ = θ,
        ξ = ξ,
    )

    ψs = collect(range(ψ_rng..., length=1001))
    ps = (1 .+ abs.(ψs)) .* (1 .- (ψs .- ψ0 ) / (coe.μ .* (1 .+ coe.θ))) ./ ( 1 .+ coe.ν * ξ / (coe.μ .* (1 .+ coe.θ)) * ( 1 .+ abs.(ψs)) )

    d_dydτ_dy = - ( 1 .+ abs.(ψs) ) + coe.μ .* (1 .+ coe.θ) * sign.(ψs) .* (1 .- (ψs .+ coe.ν * ps * ξ) / (coe.μ .* (1 .+ coe.θ)))
    stable = d_dydτ_dy .< 0

    vals, rngs = detectRanges(stable)
    
    for (k, rng) in enumerate(rngs)

        args = Dict()
        if k==1
            args[:label] = format("\$ \\xi = {:.1f} \$", -ξ)
        else
            args[:label] = nothing
        end

        args[:linestyle] = (vals[k] == 0) ? "dashed" : "solid"
        args[:color] = ["black", "red", "blue", "green", "orange", "purple"][cnt]
        ax.plot(ps[rng], ψs[rng]; args...)

    end

    global cnt += 1

end

ax.set_xlim([0, 3])
ax.set_ylim([-1, 5])
ax.legend(loc="center right", )
fig.savefig("figures/figure-stommel_bifurcation_analytical_for_AOFD_2.png", dpi=300)

plt.show()
