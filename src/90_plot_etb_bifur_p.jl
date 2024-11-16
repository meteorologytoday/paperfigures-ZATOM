using Formatting
using Roots

include("detect_ranges.jl")


ξ_ϵs = [
    (  1.0,   0.0),
    (  0.749, 0.0),
    (  0.5,   0.0),
    (  0.0,   0.0),
    (- 5.0,   0.0),
    (  0.0,   0.5),
]

ψ0 = 0.0
ψ_rng = [-4, 6]


label_positions = [
    (1.25, 3.7),
    (2.00, 3.22),
    (2.10, 2.6),
    (0.70, 1.7),
    (0.80, 0.5),
    (2.10, 1.5),
]

digits = Dict(
    :ξ => [ 0, 3, 1, 0, 0, 0],
    :ϵ => [ 0, 0, 0, 0, 0, 1],
)



μ = 3.0
ν = 1.0

println("Loading PyPlot")
println("Setting backend as Agg...")
ENV["MPLBACKEND"]="agg"
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel("\$p\$", fontsize=25)
ax.set_ylabel("\$\\Psi\$", fontsize=25)
ax.grid()

for (i, (ξ, ϵ)) in enumerate(ξ_ϵs)

    coe = (
        μ = μ,
        ν = ν,
        ξ = ξ,
    )

    ψs = collect(range(ψ_rng..., length=1001))
    ϵ_switch = (sign.(ψs) .+ 1) / 2
    ps = (1 .+ abs.(ψs)) .* (1 .- ψs / μ) ./ ( 1 .- ϵ_switch .* ϵ / μ .* ψs .- ν * ξ / μ * ( 1 .+ abs.(ψs)) )

    ys = ps ./ (1 .+ abs.(ψs))
    d_dydτ_dy = - ( 1 .+ abs.(ψs) ) + ys .* sign.(ψs) .* (μ .- ϵ .* ψs) ./ (1 .- ϵ .* ys)

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
    
    
    fmt_ξ = (digits[:ξ][i] == 0) ? ":d" : ":.$(digits[:ξ][i])f"
    fmt_ϵ = (digits[:ϵ][i] == 0) ? ":d" : ":.$(digits[:ϵ][i])f"

    ax.text(label_positions[i]..., format("\$(\\tilde{{\\xi}}, \\epsilon) = ({$fmt_ξ}, {$fmt_ϵ}) \$", ξ, ϵ), va="center", ha="center", color=color, size=13)

end

# B_1
B1 = ( 0, μ )
ax.scatter(B1..., s=20, marker="o", color="black", zorder=99)
#ax.annotate("\$ B_1 = \\left( 0, \\mu \\right) \$", xy=B1,  xycoords="data",
ax.annotate("\$ B_1 \$", xy=B1,  xycoords="data",
            xytext=(0.25, 3.5), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=13,
)

# B_2
ξ = 0.0
ϵ = 0.0
# Look numerically for the bifurcation point
α = ν * ξ / μ
β = ϵ / μ
f(ψ) = ( α + β ) * ψ^2 - 2 * (1 - α) * ψ + (μ - 1 + α + μ * β)
dfdψ(ψ) = 2 * ( α + β ) * ψ - 2 * (1 - α)
p(ψ) = (1 + abs(ψ)) * (1 - ψ / μ ) / ( 1 - ν * ξ / μ  * ( 1 + abs(ψ)) - ϵ * ψ / μ )
ψ_bifur = find_zero((f, dfdψ), 0.0, Roots.Newton())

B2 = ( p(ψ_bifur), ψ_bifur )
ax.scatter(B2..., s=20, marker="o", color="black", zorder=99)
ax.annotate("\$ B_2 \$", xy=B2,  xycoords="data",
            xytext=(1.0, 1.0), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=13,
)

# B_3
B3 = ( p(0.0), 0 )
ax.scatter(B3..., s=20, marker="o", color="black", zorder=99)
ax.annotate("\$ B_3 \$", xy=B3,  xycoords="data",
            xytext=(1.0, -0.5), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="top", fontsize=13,
)


ax.set_xlim([-0.1, 3])
ax.set_ylim([-1, 5])
ax.text(0.05, 0.95, "(a)", size=25, va="top", ha="left", transform=ax.transAxes)
#ax.legend(loc="center right",)
fig.savefig("figures/figure-etb_bifur_p.svg", dpi=300)

#plt.show()
