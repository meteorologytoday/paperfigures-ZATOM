using Formatting
using Roots
using PyCall

shp_geo = pyimport("shapely.geometry")
shp_ops = pyimport("shapely.ops")
plt_patches = pyimport("matplotlib.patches")

include("ZATOM_regimes.jl")

ξs = collect(range( -10,  0,  length=2001))
ps = collect(range( -0.1, 0.5,  length=2001))

μ = 16.1
ν = 162.0

ϕn = 70.0
ϕs = 10.0
H = 4500.0
a = 6.4e6
α_T = 2e-3
α_S = 7e-3
V   = π * H / 9 * a^2 * (sin(deg2rad(ϕn)) - sin(deg2rad(ϕs))) / 4
A_w =  5.0
A_e = 15.0

geometry_factor = 0.25 * (A_w + A_e)^2.0 / A_w / A_e 

δT_star = 25.0
S0 = 35.0

L = a * deg2rad(ϕn - ϕs)
Kh = 4e4

t_d = L^2 / Kh
t_a = V / t_d 
#t_d = 86400 * 365 * 1e3

factor_γ2p = (2 * geometry_factor * α_S * t_d * S0) / (α_T * δT_star * V)
factor_p2γ = factor_γ2p^(-1)
factor_p2γSv = factor_p2γ / 1e6

ξ_corner = -1.5
μ = 1 / (0.05e6 * factor_γ2p)
ν = - μ * (μ - 1) / ξ_corner

println("factor_γ2p = $factor_γ2p")
println("t_d = $t_d")
println("V = $V")
println("μ = $μ")
println("ν = $ν")


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

        α = - coe.ν * coe.ξ / coe.μ
     
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

ξ_idx = isfinite.(p_left_bnd .* p_right_bnd)
ξs = ξs[ξ_idx]
p_left_bnd = p_left_bnd[ξ_idx]
p_right_bnd = p_right_bnd[ξ_idx]

p_idx = isfinite.(ξ_left_bnd .* ξ_right_bnd)
ps = ps[ξ_idx]
ξ_left_bnd = ξ_left_bnd[ξ_idx]
ξ_right_bnd = ξ_right_bnd[ξ_idx]


ξ_upper_bound = μ/ν
ξ_lower_bound = - μ/ν * ( μ - 1 )

ξ_critical = μ / ( ν * ( μ+1 ) )

p_lower_bound = 1 / μ

γ_left_bnd = p_left_bnd * factor_p2γSv
γ_right_bnd = p_right_bnd * factor_p2γSv
γ_lower_bound = p_lower_bound * factor_p2γSv


println("γ_lower_bound = $(γ_lower_bound / 1e6) Sv")

fixed_ξ = zeros(Float64, length(ξs), 3)
fixed_γ = zeros(Float64, length(ps), 3)

fixed_ξ[:, 1] = ξs
fixed_ξ[:, 2] = γ_left_bnd
fixed_ξ[:, 3] = γ_right_bnd

fixed_γ[:, 1] = ps * factor_p2γSv
fixed_γ[:, 2] = ξ_left_bnd
fixed_γ[:, 3] = ξ_right_bnd

regimes["extended two-box"] = Dict(
    "label"     => "Extended\ntwo-box model",
    "label_pos" => (0.11, -0.5),
    "fixed_ξ" => fixed_ξ,
    "fixed_γ" => fixed_γ, 
)

regimes["standard"]["label_pos"] = (0.11, -1.1)
regimes["standard"]["label"] = "ZATOM"




# Creating polygons
function createPoly(ax_direction, ax_vals, l_bnds, u_bnds)

    local pts = []

    for i = 1:length(ax_vals)

        pt = [ax_vals[i], l_bnds[i]]
        if ax_direction == "vertical"
            pt[1], pt[2] = pt[2], pt[1]
        end
        push!(pts, pt)
    end

    for i = length(ax_vals):-1:1
        pt = [ax_vals[i], u_bnds[i]]
        if ax_direction == "vertical"
            pt[1], pt[2] = pt[2], pt[1]
        end
        push!(pts, pt)
    end

    return shp_geo.Polygon(pts)
end

function shpPoly2MatplotPoly(poly, kwarg)
    poly_x , poly_y = poly.exterior.xy
    poly_xy = zeros(Float64, length(poly_x), 2)
    poly_xy[:, 1] = poly_x
    poly_xy[:, 2] = poly_y
    return plt_patches.Polygon(poly_xy; kwarg...)
end

println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fcs = [ "dodgerblue", "orangered"]
fcs = [ "none", "none"]
ecs = [ "blue", "red"]
hatches = [ "..", "//" ]

fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel("\$\\gamma\$ [Sv]", fontsize=25)
ax.set_ylabel("\$\\xi\$", fontsize=25)
ax.grid()
ax.set_ylim([-1.6, -0.4])
ax.set_xlim([0.04, 0.15])

#ax.fill_betweenx(ξs, γ_left_bnd, γ_right_bnd, facecolor="blue", edgecolor="blue",       hatch="..", alpha=0.8, linewidth=1, zorder=10)#, label="Folding along fixed \$\\xi\$")
#ax.fill_between(ps * factor_p2γSv, ξ_left_bnd, ξ_right_bnd, facecolor="none",  edgecolor="orangered",  hatch="//", alpha=0.8, linewidth=1, zorder=10)#, label="Folding along fixed \$p\$")

for (k, key) in enumerate(["standard", "extended two-box"])

    println("Plotting regime = $key")

    data = regimes[key]

    fixed_ξ = data["fixed_ξ"]
    fixed_γ = data["fixed_γ"]

    label = data["label"]

    # creating polygons
    poly_ξ = createPoly("vertical",   fixed_ξ[:, 1], fixed_ξ[:, 2], fixed_ξ[:, 3])
    poly_γ = createPoly("horizontal", fixed_γ[:, 1], fixed_γ[:, 2], fixed_γ[:, 3])
    merged_poly = shp_ops.cascaded_union([poly_ξ, poly_γ])
    merged_poly = shpPoly2MatplotPoly(merged_poly, Dict(
        :ec     => ecs[k],
        :fc     => fcs[k],
        :alpha  => 1.0,
        :zorder => 20,
        :hatch  => hatches[k],
    ))


    #ax.fill_betweenx(fixed_ξ[:, 1], fixed_ξ[:, 2], fixed_ξ[:, 3], hatch="..", facecolor="none",  edgecolor=colors["fixed_ξ"], alpha=0.8, linewidth=1, zorder=10, label="[$label] Folding along fixed \$\\xi\$")
    #ax.fill_between(fixed_γ[:, 1], fixed_γ[:, 2], fixed_γ[:, 3], hatch="//",  facecolor="none",  edgecolor=colors["fixed_γ"], alpha=0.8, linewidth=1, zorder=10, label="[$label] Folding along fixed \$\\gamma\$")

    ax.add_patch(merged_poly)

    ax.text(data["label_pos"]..., data["label"], size=12, ha="center", va="center", color=ecs[k])
end



ax.plot(ax.get_xlim(), [ξ_upper_bound, ξ_upper_bound], ls="dashed", color="#aaaaaa")
ax.plot(ax.get_xlim(), [ξ_lower_bound, ξ_lower_bound], ls="dashed", color="#aaaaaa")
ax.plot([γ_lower_bound, γ_lower_bound], ax.get_ylim(), ls="dashed", color="#aaaaaa")

#ax.set_title("(b)")
#ax.legend()

fig.savefig("figures/regime_diagrams_comparison.png", dpi=300)

plt.show()
