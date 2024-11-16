using Formatting
include("detect_ranges.jl")


ps = [ 0.1, 0.5, 1.0, 2.0, 1e5, 1.0]
ϵs = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.5]

digits = Dict(
    :p => [1, 1, 0, 0, 0, 0],
    :ϵ => [0, 0, 0, 0, 0, 1],
)


ψ0 = 0.0
ψ_rng = [-4, 6]

μ = 3.0
ν = 1.0

label_positions = [
    (-4.50, 2.75),
    (-4.50, 0.4),
    (-1.80, 0.4),
    ( 0.75, 0.4),
    ( 2.50, 1.5),
    (-2.00, 1.3),
]


println("Loading PyPlot")
println("Setting backend as Agg...")
ENV["MPLBACKEND"]="agg"
using PyPlot
plt = PyPlot
println("Done")


fig, ax = plt.subplots(1, 1, constrained_layout=true)



ax.set_xlabel("\$\\tilde{\\xi}\$", fontsize=25)
ax.set_ylabel("\$\\Psi\$", fontsize=25)
ax.grid()

for (i, p) in enumerate(ps)

    coe = (
        c = 3.6e3,
        μ = μ,
        ν = ν,
        p = p,
        ϵ = ϵs[i],
    )


    ψs = collect(range(ψ_rng..., length=1001))
    
    ηs = (sign.(ψs) .+ 1) / 2
    ξs = coe.μ / coe.ν * ( (1 .- ηs .* coe.ϵ .* ψs ./ coe.μ ) ./ (1 .+ abs.(ψs)) .- (1 .- ψs / coe.μ) / coe.p )

    ys = coe.p ./ (1 .+ abs.(ψs))
    d_dydτ_dy = - ( 1 .+ abs.(ψs) ) + ys .* sign.(ψs) .* (coe.μ .- coe.ϵ .* ψs) ./ (1 .- coe.ϵ .* ys)
    
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

        fmt_p = (digits[:p][i] == 0) ? ":d" : ":.$(digits[:p][i])f"
        fmt_ϵ = (digits[:ϵ][i] == 0) ? ":d" : ":.$(digits[:ϵ][i])f"

        if p < 1e3
            fmt 
            txt = format("\$(p, \\epsilon) = ({$fmt_p}, {$fmt_ϵ}) \$", p, coe.ϵ)
        else
            txt = format("\$(p, \\epsilon) = (\\infty, {$fmt_ϵ}) \$", coe.ϵ)
        end
        ax.text(label_positions[i]..., txt, va="center", ha="center", color=color, size=13)
    end

end

p = 1.0
ξ(ψ) = μ / ν * ( 1 / (1 + abs(ψ)) - (1 - ψ / μ) / p )

# C1
# ξ, ψ
C1 = ( μ / (1 + μ) / ν , μ )
ax.scatter(C1..., s=20, marker="o", color="black", zorder=99)
#ax.annotate("\$ C_1 = \\left( \\frac{\\mu}{\\left( 1 + \\mu \\right) \\nu }, \\mu \\right) \$", xy=C1,  xycoords="data",
ax.annotate("\$ C_1 \$", xy=C1,  xycoords="data",
            xytext=(0.3, 3.5), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=13,
)

# C2
ψ_bifur = (p * μ)^0.5 - 1
C2 = ( ξ(ψ_bifur), ψ_bifur )
ax.scatter(C2..., s=20, marker="o", color="black", zorder=99)
ax.annotate("\$ C_2 \$", xy=C2,  xycoords="data",
            xytext=(-1.5, ψ_bifur), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=13,
)

# C3
C3 = ( μ / ν * (1 - 1 / p), 0 )
ax.scatter(C3..., s=20, marker="o", color="black", zorder=99)
#ax.annotate("\$ C_3 = \\left( \\frac{\\mu}{\\nu} \\left( 1 - \\frac{1}{p} \\right), 0 \\right) \$", xy=C3,  xycoords="data",
ax.annotate("\$ C_3 \$", xy=C3,  xycoords="data",
            xytext=(-1.5, 0), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=13,
)


#=
# Plot the line connecting C2
ξ(ψ, p) = μ / ν * ( 1 / (1 + abs(ψ)) - (1 - ψ / μ) / p )
ps = range(0, 100, length=1001) |> collect
ψ_bifurs = (ps .* μ).^0.5 .- 1
ξs = [ξ(ψ_bifurs[i], ps[i]) for i = 1:length(ps)]

ax.plot(ξs, ψ_bifurs, "--", color="#bbbbbb", zorder=99)
#ax.annotate("\$ C_3 = \\left( \\frac{\\mu}{\\nu} \\left( 1 - \\frac{1}{p} \\right), 0 \\right) \$", xy=C3,  xycoords="data",
ax.annotate("\$ C_3 \$", xy=C3,  xycoords="data",
            xytext=(-1.5, 0), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=13,
)
=#

ax.set_xlim([-6, 4])
ax.set_ylim([-1, 5])
ax.text(0.05, 0.95, "(b)", size=25, va="top", ha="left", transform=ax.transAxes)
#ax.legend(loc="center right",)

fig.savefig("figures/figure-etb_bifur_xi.svg", dpi=300)

#plt.show()
