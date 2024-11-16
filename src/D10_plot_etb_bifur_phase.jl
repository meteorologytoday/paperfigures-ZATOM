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

       
        α = ν * ξ / μ
        β = ϵ / μ
        f(ψ) = ( α + β ) * ψ^2 - 2 * (1 - α) * ψ + (μ - 1 + α + μ * β)
        dfdψ(ψ) = 2 * ( α + β ) * ψ - 2 * (1 - α)
        p(ψ) = ( 1 + abs(ψ) ) * ( 1 - ψ / μ ) / ( 1 - ν * ξ / μ  * ( 1 + abs(ψ) ) - ϵ * ψ / μ )

        Δ = (1 - α)^2 - (α+β) * ((μ-1)*(1-α)+μ*(α+β))

        #Δ = (1+α)^2 - α * (1 + α - coe.μ )
        #Δ2 = coe.μ ^2 - coe.μ + coe.ν * ξ
        
        local ψ_left, ψ_right, p_left, p_right

        if Δ >= 0 

            ψ_left = 0.0
            ψ_right = find_zero((f, dfdψ), 0.0, Roots.Newton()) 

            if ψ_right < 0  # it is positive definite
                continue
            end

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

    ps[(ps .* ( μ + ϵ ) ) .< 1] .= NaN

    Ψ_c  = -1 .+ (ps .* (μ + ϵ) ).^0.5   # Critical Ψ
    #- (1 .+ ps .* ϵ ./ 2) .+ ( (1 .+ ps .* ϵ ./ 2).^2 .+ ps .* μ .- 1 ).^0.5   # Critical Ψ

    ξ_left_bnd  = μ / ν * ( (1 .- ϵ .* Ψ_c ./ μ) ./ (1 .+ abs.(Ψ_c)) - (1 .- Ψ_c / μ) ./ ps )
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
ξ_lower_bound = - μ/ν * ( μ + ϵ - 1 )

ξ_critical = μ * (1 - ϵ) / ( ν * ( μ + 1 ) )

p_lower_bound = 1 / ( μ + ϵ )


println("Loading PyPlot")
println("Setting backend as Agg...")
ENV["MPLBACKEND"]="agg"
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel("\$p\$", fontsize=25)
ax.set_ylabel("\$\\tilde{\\xi}\$", fontsize=25)
ax.grid(alpha=0.5)
ax.set_ylim([-8, 4.0])
ax.set_xlim([0.0, 5])

# Just show one direction is sufficient
#ax.fill_betweenx(ξs, p_left_bnd, p_right_bnd, facecolor="none", edgecolor="blue",       hatch="..", alpha=0.8, linewidth=1, zorder=10, label="Folding along fixed \$\\xi\$")
ax.fill_between(ps, ξ_left_bnd, ξ_right_bnd, facecolor="orangered",  edgecolor="red",  alpha=0.5, linewidth=1, zorder=10, label="Folding along fixed \$p\$")

ax.plot(ax.get_xlim(), [ξ_upper_bound, ξ_upper_bound], ls="dashed", color="#000000")

ax.text(1.5, ξ_upper_bound + 0.2, "\$\\tilde{\\xi} = \\tilde{\\xi}_u = \\frac{\\mu}{\\nu} \$", ha="center", va="bottom", color="#000000")


#ax.annotate("\$ \\tilde{\\xi} = \\frac{\\mu \\left(1 - \\epsilon \\right)}{\\nu \\left( \\mu + 1 \\right)} \$", xy=(1.1, ξ_critical),  xycoords="data",
#            xytext=(0.5, ξ_critical), textcoords="data",
#            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
#            horizontalalignment="center", verticalalignment="center", fontsize=13,
#)


ax.scatter(p_lower_bound, ξ_lower_bound, s=20, color="black", zorder=99)
ax.annotate("\$ P = \\left(p_c, \\tilde{\\xi}_c \\right) = \\left( \\frac{1}{\\mu + \\epsilon}, - \\frac{\\mu  \\left( \\mu + \\epsilon - 1 \\right) }{\\nu } \\right) \$",
            xy=(p_lower_bound, ξ_lower_bound),  xycoords="data",
            xytext=(p_lower_bound + 0.5, ξ_lower_bound), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "headlength" => 5.0, "width" => 0.5),
            horizontalalignment="left", verticalalignment="center", fontsize=13,
)


#ax.set_title("(b)")

#ax.legend(loc="center right")

fig.savefig("figures/figure-etb_bifur_phase.svg", dpi=300)

println("Showing figure...")
#plt.show()
