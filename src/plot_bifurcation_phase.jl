using Formatting
include("detect_ranges.jl")


ξs = collect(range(-1, 5,  length=51))

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
    )

    α = coe.ν * coe.ξ / coe.μ
    existence_criteria_satisfied = (1 + α - coe.μ) / α < 0.0
    p(ψ) = (1 + abs(ψ)) * (1 - ψ / coe.μ) / ( 1 + coe.ν * coe.ξ / coe.μ * ( 1 + abs(ψ)) )

    local ψ_left, ψ_right, p_left, p_right

    if existence_criteria_satisfied 

        ψ_left = 0.0
        α_tmp = (1 + α) / α
        ψ_right = - α_tmp + (α_tmp^2 + coe.μ / α - α_tmp )^0.5

        p_left  = p(ψ_left)
        p_right = p(ψ_right)

        if p_right < p_left
            println("p_left = $p_left , p_right = $p_right")
            throw(ErrorException("We have p_left > p_right which is wrong."))
        end
   
        p_left_bnd[j]  = p_left 
        p_right_bnd[j] = p_right

    end

end

println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel("\$p\$", fontsize=25)
ax.set_ylabel("\$\\xi\$", fontsize=25)
ax.grid()
ax.invert_yaxis()

ax.plot(p_left_bnd, ξs, color="red")
ax.plot(p_right_bnd, ξs, color="blue")

plt.show()
