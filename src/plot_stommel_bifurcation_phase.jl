using Formatting
using Roots


ξs = collect(range(-10, 5,  length=2001))
ps = collect(range( -1, 5,  length=2001))
μ = 3.0
ν = 1.0
ϵ = 0.5

function generate_bnds_p(ξs)

    p_left_bnd  = zeros(Float64, length(ξs))
    p_right_bnd = zeros(Float64, length(ξs))

    p_left_bnd .= NaN
    p_right_bnd .= NaN

    for (j, ξ) in enumerate(ξs)

        coe = (
            μ = μ,
            ν = ν,
            ξ = ξ,
        )

        α = coe.ν * coe.ξ / coe.μ
     
        Δ = (1+α)^2 - α * (1 + α - coe.μ )
        Δ2 = coe.μ ^2 - coe.μ + coe.ν * ξ
        
        f(ψ) = α * ψ^2 + 2 * (1 + α) * ψ + (1 + α - coe.μ)
        dfdψ(ψ) = 2 * α * ψ + 2 * (1 + α)

        p(ψ) = (1 + abs(ψ)) * (1 - ψ / coe.μ ) / ( 1 - coe.ν * ξ / coe.μ  * ( 1 + abs(ψ)) )

        local ψ_left, ψ_right, p_left, p_right

        if Δ >= 0 && Δ2 > 0

            ψ_left = 0.0
            α_tmp = (1 + α) / α
            ψ_right = find_zero((f, dfdψ), 0.0, Roots.Newton()) 



            p_left  = p(ψ_left)
            p_right = p(ψ_right)
            
            #println("ξ = $ξ , α = $α , Δ = $Δ , ψ_right = $ψ_right , p_right = $p_right")

            if p_right < p_left
                println("p_left = $p_left , p_right = $p_right")
                throw(ErrorException("We have p_left > p_right which is wrong."))
            end
       
            p_left_bnd[j]  = p_left 
            p_right_bnd[j] = p_right

        end

    end

    return p_left_bnd, p_right_bnd

end

function generate_bnds_ξ(ps)

    ps = copy(ps)

    ps[(ps .* μ) .< 1] .= NaN

    Ψ_c  = ( ps .* μ ).^0.5 .- 1   # Critical Ψ

    ξ_left_bnd  = μ / ν * ( 1 ./ (1 .+ abs.(Ψ_c)) - (1 .- Ψ_c / μ) ./ ps )
    ξ_right_bnd = μ / ν .* ( 1 .- 1 ./ ps )

    no_solution_idx = ξ_left_bnd .> ξ_right_bnd

    ξ_left_bnd[no_solution_idx] .= NaN
    ξ_right_bnd[no_solution_idx] .= NaN

    return ξ_left_bnd, ξ_right_bnd

end

println("Compute boundaries...")
p_left_bnd, p_right_bnd = generate_bnds_p(ξs)
ξ_left_bnd, ξ_right_bnd = generate_bnds_ξ(ps)

ξ_upper_bound = μ/ν
ξ_lower_bound = - μ/ν * ( μ - 1 )

ξ_critical = μ / ( ν * ( μ+1 ) )

p_lower_bound = 1 / μ


println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel("\$p\$", fontsize=25)
ax.set_ylabel("\$\\xi\$", fontsize=25)
ax.grid()
ax.set_ylim([-7, 4.0])
ax.set_xlim([0.0, 5])


ax.fill_betweenx(ξs, p_left_bnd, p_right_bnd, facecolor="none", edgecolor="blue",       hatch="..", alpha=0.8, linewidth=1, zorder=10, label="Folding along fixed \$\\xi\$")
ax.fill_between(ps, ξ_left_bnd, ξ_right_bnd, facecolor="none",  edgecolor="orangered",  hatch="//", alpha=0.8, linewidth=1, zorder=10, label="Folding along fixed \$p\$")

ax.plot(ax.get_xlim(), [ξ_upper_bound, ξ_upper_bound], ls="dashed", color="#aaaaaa")
ax.plot(ax.get_xlim(), [ξ_lower_bound, ξ_lower_bound], ls="dashed", color="#aaaaaa")

ax.text(1.5, ξ_upper_bound + 0.2, "\$\\xi = \\frac{\\mu}{\\nu} \$", ha="center", va="bottom", color="#555555")
ax.text(1.5, ξ_lower_bound + 0.2, "\$\\xi = - \\frac{\\mu  \\left( \\mu - 1 \\right) }{\\nu }\$", ha="center", va="bottom", color="#555555")

ax.plot([p_lower_bound, p_lower_bound], ax.get_ylim(), ls="dashed", color="#aaaaaa")
ax.text(p_lower_bound + 0.1, 1.8, "\$ p = \\frac{1}{\\mu} \$", ha="left", va="center", rotation=0, color="#555555")


ax.annotate("\$ \\xi = \\frac{\\mu}{\\nu \\left( \\mu + 1 \\right)} \$", xy=(1.2, ξ_critical),  xycoords="data",
            xytext=(0.5, ξ_critical), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=13,
)

ax.scatter(p_lower_bound, ξ_lower_bound, s=20, color="black", zorder=99)
ax.annotate("\$ P \$", xy=(p_lower_bound, ξ_lower_bound),  xycoords="data",
            xytext=(p_lower_bound + 0.5, ξ_lower_bound), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=13,
)


#ax.set_title("(b)")

ax.legend()

fig.savefig("figures/figure-stommel_bifurcation_phase.png", dpi=300)

plt.show()
