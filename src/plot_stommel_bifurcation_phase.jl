using Formatting
using Roots


ξs = collect(range(-5, 10,  length=1001))


function generate_bnds(ξs, θ)

    println("Doing θ = $θ")

    p_left_bnd  = zeros(Float64, length(ξs))
    p_right_bnd = zeros(Float64, length(ξs))

    p_left_bnd .= NaN
    p_right_bnd .= NaN

    for (j, ξ) in enumerate(ξs)

        coe = (
            c = 3.6e3,
            μ = 3.0,
            ν = 1.0,
            ξ = ξ,
            θ = θ,
        )

        α = coe.ν * coe.ξ / ( coe.μ * (1 + coe.θ) )
     
        Δ = (1+α)^2 - α * (1 + α - coe.μ * (1 + coe.θ) )
        Δ2 = ( coe.μ * (1 + coe.θ) )^2 - coe.μ * (1 + coe.θ) - coe.ν * ξ
        
        f(ψ) = α * ψ^2 + 2 * (1 + α) * ψ + (1 + α - coe.μ * (1 + coe.θ))
        dfdψ(ψ) = 2 * α * ψ + 2 * (1 + α)

        p(ψ) = (1 + abs(ψ)) * (1 - ψ / ( coe.μ * (1 + coe.θ) )) / ( 1 + coe.ν * ξ / ( coe.μ * (1 + coe.θ) ) * ( 1 + abs(ψ)) )

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

println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel("\$p\$", fontsize=25)
ax.set_ylabel("\$\\xi\$", fontsize=25)
ax.grid()
ax.set_ylim([-1, 6.2])
ax.set_xlim([0.0, 6])
ax.invert_yaxis()


#ax.plot(p_left_bnd, ξs, color="red")
#ax.plot(p_right_bnd, ξs, color="blue")
for (i, θ) in enumerate([ 0.5, 0 ])

    facecolor = ["dodgerblue", "lightsalmon", "dodgerblue", "lightgreen"][i]
    edgecolor = ["blue", "orangered", "blue", "green"][i]
    

    p_left_bnd, p_right_bnd = generate_bnds(ξs, θ)
    ax.fill_betweenx(ξs, p_left_bnd, p_right_bnd, facecolor=facecolor, edgecolor=edgecolor, alpha=0.8, linewidth=1, zorder=10, label="\$\\theta=$θ\$")
end
ax.legend(loc="lower right")
ax.set_title("Multiple equilibria phase diagram")
fig.savefig("figures/figure-stommel_bifurcation_phase.png", dpi=300)

plt.show()
