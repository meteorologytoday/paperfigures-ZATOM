using Formatting
using Roots
using PyCall

shp_geo = pyimport("shapely.geometry")
shp_ops = pyimport("shapely.ops")
plt_patches = pyimport("matplotlib.patches")

include("ZATOM_regimes.jl")

mutable struct ETBDims
    
    t_d :: Float64
    t_R :: Float64
    V   :: Float64
    δT_star :: Float64
    G   :: Float64
    μ   :: Float64
    ν   :: Float64

    factors :: Any

    function ETBDims(t_d, t_R, V, δT_star, G, ξ_corner, γ_corner)

        factor_γ2p = (2 * (1 + G) * α_S * t_d * S0) / (α_T * δT_star * V)
        factor_p2γ = factor_γ2p^(-1)
        factor_p2γSv = factor_p2γ / 1e6
        μ = 1 / (γ_corner * factor_γ2p)
        ν = - μ * (μ - 1) / ξ_corner

        println("t_d = $(t_d / 86400/365) yr.")
        println("factor_γ2p = $factor_γ2p")
        println("μ = $μ")
        println("ν = $ν")


        return new(

            t_d, t_R, V, δT_star, G, μ, ν,
 
            (
                γ2p   = factor_γ2p,
                p2γ   = factor_p2γ,
                p2γSv = factor_p2γSv,
            ),

        )

    end


end


ξs_all = collect(range( -5,  0.5,  length=10001))
ps_all = collect(range( -0.1, 0.5,  length=10001))

ϕn = 70.0
ϕs = 10.0
H = 4500.0
a = 6.4e6
α_T = 2e-3
α_S = 7e-3
V   = π * H / 9 * a^2 * (sin(deg2rad(ϕn)) - sin(deg2rad(ϕs))) / 2
A_w =  5.0
A_e = 15.0

G = 0.5 * (A_w / (2 * A_e) + A_e / (2 * A_w) - 1)

δT_star = 25.0
S0 = 35.0

L = a * deg2rad(ϕn - ϕs)
Kh = 4e4

t_d = L^2 / Kh
t_a = V / t_d 
#t_d = 86400 * 365 * 1e3
t_R = 10 * 86400.0

shifted_ξ = parse(Float64, ARGS[1])
shifted_γ = parse(Float64, ARGS[2])
ξ_corner = -1.5   + shifted_ξ
γ_corner = (0.05 + shifted_γ) * 1e6

if ARGS[3] == "plot_shifted"
    plot_shifted = true
else
    plot_shifted = false
end

println("L = $L")
println("t_d = $t_d")
println("V = $V")

etb_dims = Dict(

    "etb_zatom" => ETBDims(
        t_d,
        t_R,
        V,
        δT_star,
        G,
        ξ_corner,
        γ_corner,
    ),


    "etb_short_t_d" => ETBDims(
        1 * 86400 * 365.0,
        t_R,
        V,
        δT_star,
        G,
        ξ_corner,
        γ_corner,
    ),

    "etb_longer_t_d" => ETBDims(
        350 * 86400 * 365.0,
        t_R,
        V,
        δT_star,
        G,
        ξ_corner,
        γ_corner,
    ),

)


function generate_bnds_p(ξs, ed :: ETBDims)

    p_left_bnd  = zeros(Float64, length(ξs))
    p_right_bnd = zeros(Float64, length(ξs))

    p_left_bnd .= NaN
    p_right_bnd .= NaN

    for (j, ξ) in enumerate(ξs)

        coe = (
            μ = ed.μ,
            ν = ed.ν,
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

    ξ_idx = isfinite.(p_left_bnd .* p_right_bnd)
    ξs = ξs[ξ_idx]
    p_left_bnd = p_left_bnd[ξ_idx]
    p_right_bnd = p_right_bnd[ξ_idx]


    return ξs, p_left_bnd, p_right_bnd

end

function generate_bnds_ξ(ps, ed :: ETBDims)

    ps = copy(ps)

    ps[(ps .* ed.μ) .< 1] .= NaN

    Ψ_c  = ( ps .* ed.μ ).^0.5 .- 1   # Critical Ψ

    ξ_left_bnd  = ed.μ / ed.ν * ( 1 ./ (1 .+ abs.(Ψ_c)) - (1 .- Ψ_c / ed.μ) ./ ps )
    ξ_right_bnd = ed.μ / ed.ν .* ( 1 .- 1 ./ ps )

    no_solution_idx = ξ_left_bnd .> ξ_right_bnd

    ξ_left_bnd[no_solution_idx] .= NaN
    ξ_right_bnd[no_solution_idx] .= NaN

    p_idx = isfinite.(ξ_left_bnd .* ξ_right_bnd)
    ps = ps[p_idx]
    ξ_left_bnd = ξ_left_bnd[p_idx]
    ξ_right_bnd = ξ_right_bnd[p_idx]

    return ps, ξ_left_bnd, ξ_right_bnd

end

println("Compute boundaries...")

data = Dict()

for (k, ed) in etb_dims 

    ξs, p_left_bnd, p_right_bnd = generate_bnds_p(ξs_all, ed)
    ps, ξ_left_bnd, ξ_right_bnd = generate_bnds_ξ(ps_all, ed)

    println("Length of ξs = ", length(ξs))
    println("Length of ps = ", length(ps))
    
    ξ_upper_bound = ed.μ / ed.ν
    ξ_lower_bound = - ed.μ / ed.ν * ( ed.μ - 1 )

    ξ_critical = ed.μ / ( ed.ν * ( ed.μ + 1 ) )

    p_lower_bound = 1 / ed.μ

    γ_left_bnd = p_left_bnd * ed.factors.p2γSv
    γ_right_bnd = p_right_bnd * ed.factors.p2γSv
    γ_lower_bound = p_lower_bound * ed.factors.p2γSv
    γs = ps * ed.factors.p2γSv

    data[k] = Dict(
        "ξ_bnds" => (ξ_left_bnd, ξ_right_bnd),
        "p_bnds" => (p_left_bnd, p_right_bnd),
        "γ_bnds" => (γ_left_bnd, γ_right_bnd),
        "ps" => ps,
        "γs" => γs,
        "ξs" => ξs,
        "γ_lower_bound" => γ_lower_bound,
    )

end


regimes["etb_zatom"] = Dict(
        "label"     => "ETBM \$t_d=35\\mathrm{yr}\$",
        "label_pos" => (0.134, -0.4),
)

regimes["etb_short_t_d"] = Dict(
        "label"     => "ETBM \$t_d=1\\mathrm{yr}\$",
        "label_pos" => (0.16, -0.6),
)


regimes["etb_longer_t_d"] = Dict(
        "label"     => "ETBM \$t_d=350\\mathrm{yr}\$",
        "label_pos" => (0.12, -0.1),
)

for k in keys(data)

    ξs = data[k]["ξs"]
    γs = data[k]["γs"]

    fixed_ξ = zeros(Float64, length(ξs), 3)
    fixed_γ = zeros(Float64, length(γs), 3)

    println(k)

    fixed_ξ[:, 1] = ξs

    fixed_ξ[:, 2] = data[k]["γ_bnds"][1]
    fixed_ξ[:, 3] = data[k]["γ_bnds"][2]

    fixed_γ[:, 1] = γs
    fixed_γ[:, 2] = data[k]["ξ_bnds"][1]
    fixed_γ[:, 3] = data[k]["ξ_bnds"][2]

    merge!(regimes[k], Dict(
        "fixed_ξ" => fixed_ξ,
        "fixed_γ" => fixed_γ, 
    ))
    

end

regimes["standard"]["label_pos"] = (0.11, -1.1)
regimes["standard"]["label"] = "ZATOM"

regimes["standard_adjusted"] = deepcopy(regimes["standard"])
regimes["standard_adjusted"]["label"] = "ZATOM adjusted"
regimes["standard_adjusted"]["fixed_ξ"][:, 1]     .+= shifted_ξ
regimes["standard_adjusted"]["fixed_ξ"][:, 2:3]   .+= shifted_γ
regimes["standard_adjusted"]["fixed_γ"][:, 1]     .+= shifted_γ
regimes["standard_adjusted"]["fixed_γ"][:, 2:3]   .+= shifted_ξ

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

plot_cases = ["etb_longer_t_d", "etb_zatom", "etb_short_t_d", "standard"]
fcs = [ "none", "none", "none", "none"][end:-1:1]
ecs = [ "blue", "red", "orange", "green" ][end:-1:1]
hatches = [ "..", "//", "//", "//" ][end:-1:1]

if plot_shifted
    push!(fcs, "none")
    push!(ecs, "#aaaaaa")
    push!(hatches, "..")
    push!(plot_cases, "standard_adjusted")
end


fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel("\$\\gamma\$ [Sv]", fontsize=25)
ax.set_ylabel("\$\\xi\$", fontsize=25)
ax.grid()
ax.set_ylim([-1.6, 0.1])
ax.set_xlim([0.04, 0.2])

#ax.fill_betweenx(ξs, γ_left_bnd, γ_right_bnd, facecolor="blue", edgecolor="blue",       hatch="..", alpha=0.8, linewidth=1, zorder=10)#, label="Folding along fixed \$\\xi\$")
#ax.fill_between(ps * factor_p2γSv, ξ_left_bnd, ξ_right_bnd, facecolor="none",  edgecolor="orangered",  hatch="//", alpha=0.8, linewidth=1, zorder=10)#, label="Folding along fixed \$p\$")

#for (k, key) in enumerate(["standard", "etb_short_t_d", "etb_zatom", "etb_longer_t_d"])
for (k, key) in enumerate(plot_cases)

    println("Plotting regime = $key")

    regime = regimes[key]

    fixed_ξ = regime["fixed_ξ"]
    fixed_γ = regime["fixed_γ"]

    label = regime["label"]

    # creating polygons
    poly_ξ = createPoly("vertical",   fixed_ξ[:, 1], fixed_ξ[:, 2], fixed_ξ[:, 3])
    poly_γ = createPoly("horizontal", fixed_γ[:, 1], fixed_γ[:, 2], fixed_γ[:, 3])
    merged_poly = shp_ops.unary_union([poly_ξ, poly_γ])
    merged_poly = shpPoly2MatplotPoly(merged_poly, Dict(
        :ec     => ecs[k],
        :fc     => fcs[k],
        :alpha  => 1.0,
        :zorder => 20,
        :hatch  => hatches[k],
        :label  => regime["label"],
    ))


    #ax.fill_betweenx(fixed_ξ[:, 1], fixed_ξ[:, 2], fixed_ξ[:, 3], hatch="..", facecolor="none",  edgecolor=colors["fixed_ξ"], alpha=0.8, linewidth=1, zorder=10, label="[$label] Folding along fixed \$\\xi\$")
    #ax.fill_between(fixed_γ[:, 1], fixed_γ[:, 2], fixed_γ[:, 3], hatch="//",  facecolor="none",  edgecolor=colors["fixed_γ"], alpha=0.8, linewidth=1, zorder=10, label="[$label] Folding along fixed \$\\gamma\$")

    ax.add_patch(merged_poly)

    if length(ARGS) >= 4
        ax.set_title(ARGS[4])
    end

    #ax.text(regime["label_pos"]..., regime["label"], size=12, ha="center", va="center", color=ecs[k])
end



#ax.plot(ax.get_xlim(), [ξ_upper_bound, ξ_upper_bound], ls="dashed", color="#aaaaaa")
#ax.plot(ax.get_xlim(), [ξ_lower_bound, ξ_lower_bound], ls="dashed", color="#aaaaaa")
#ax.plot([γ_lower_bound, γ_lower_bound], ax.get_ylim(), ls="dashed", color="#aaaaaa")

#ax.set_title("(b)")
ax.legend(loc="lower right")

fig.savefig(format("figures/regime_diagrams_comparison_xi_{:.2f}_gamma_{:.2f}.png", shifted_ξ, shifted_γ), dpi=300)

plt.show()
