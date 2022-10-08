using Formatting
using Roots
using PyCall

shp_geo = pyimport("shapely.geometry")
shp_ops = pyimport("shapely.ops")
plt_patches = pyimport("matplotlib.patches")

include("ZATOM_regimes.jl")

mutable struct ETBDims
    
    t_d :: Float64
    t_r :: Float64
    V   :: Float64
    δT_star :: Float64

    p_cusp  :: Float64
    Δξ_fold :: Float64
    Δξ_arc  :: Float64
    ξ0  :: Float64

    μ   :: Float64
    ν   :: Float64
    ϵ   :: Float64


    factors :: Any

    function ETBDims(t_d, t_r, V, δT_star, γ_cusp, Δξ_fold, Δξ_arc, ξ0)

        factor_γ2p = (2 * α_S * t_d * S0) / (α_T * δT_star * V)
        factor_p2γ = factor_γ2p^(-1)
        factor_p2γSv = factor_p2γ / 1e6
        
        p_cusp = γ_cusp * factor_γ2p

        r = Δξ_fold / Δξ_arc
        μ = ( ( 1 / p_cusp - 1) / r)^0.5
        ν = μ^2 / Δξ_arc
        ϵ = Δξ_fold * ν + 1 - μ


        println("V = $(V/1e15) * 1e15 m^3")
        println("t_d = $(t_d / 86400/365) yr.")
        println("factor_γ2p = $factor_γ2p")
        println("μ = $μ")
        println("ν = $ν")
        println("ϵ = $ϵ")

        ψ0 = μ * V / t_d
        println("ψ0 = $(ψ0 / 1e6) Sv")

        return new(

            t_d, t_r, V, δT_star, p_cusp, Δξ_fold, Δξ_arc, ξ0, μ, ν, ϵ,
 
            (
                γ2p   = factor_γ2p,
                p2γ   = factor_p2γ,
                p2γSv = factor_p2γSv,
            ),

        )

    end


end


ξs_all = collect(range( -5,  0.5,  length=1001))
ps_all = collect(range( -0.1, 5.0,  length=1001))

ϕn = 70.0
ϕs = 10.0
H = 4500.0
a = 6.4e6
α_T = 2e-3
α_S = 7e-3
V   = ((20 / 360) * 2π * H * a^2 * (sin(deg2rad(ϕn)) - sin(deg2rad(ϕs))) / 2)  *  (800/4500) * (5/20)

δT_star = 25.0
S0 = 35.0

L = a * deg2rad(ϕn - ϕs)
Kh = 4e4

t_d = L^2 / Kh 
t_a = V / t_d 
t_r = 10 * 86400.0

Q = t_d / t_r

cusp_pt = (γ=0.05e6, ξ=-1.5)

γ_cusp  = cusp_pt.γ
Δξ_fold = 0.35
Δξ_arc  = 0.9

ξ0 = -0.9

γ_corner = cusp_pt[1] * 1e6
ξ_corner = cusp_pt[2]

println("L = $L")
println("t_d = $t_d")
println("V = $V")


etb_dims = Dict(

    "etb_zatom" => ETBDims(
        t_d,
        t_r,
        V,
        δT_star,
        γ_cusp,
        Δξ_fold,
        Δξ_arc,
        ξ0,
    ),

)


function generate_bnds_p(ξs, ed :: ETBDims)

    p_left_bnd  = zeros(Float64, length(ξs))
    p_right_bnd = zeros(Float64, length(ξs))

    p_left_bnd .= NaN
    p_right_bnd .= NaN

    for (j, ξ) in enumerate(ξs)

        α = ed.ν * ξ / ed.μ
        β = ed.ϵ / ed.μ
        f(ψ) = ( α + β ) * ψ^2 - 2 * (1 - α) * ψ + (ed.μ - 1 + α + ed.μ * β)
        dfdψ(ψ) = 2 * ( α + β ) * ψ - 2 * (1 - α)
        p(ψ) = ( 1 + abs(ψ) ) * ( 1 - ψ / ed.μ ) / ( 1 - α  * ( 1 + abs(ψ) ) - ψ * β )

        Δ = (1 - α)^2 - (α + β) * ( (ed.μ - 1) * (1 - α) + ed.μ * (α + β) )

        local ψ_left, ψ_right, p_left, p_right

        if Δ >= 0

            ψ_left = 0.0
            ψ_right = find_zero((f, dfdψ), 0.0, Roots.Newton()) 

            if ψ_right < 0  # it is positive definite
                continue
            end


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

    ξ_idx = isfinite.(p_left_bnd .* p_right_bnd)
    ξs = ξs[ξ_idx]
    p_left_bnd = p_left_bnd[ξ_idx]
    p_right_bnd = p_right_bnd[ξ_idx]


    return ξs, p_left_bnd, p_right_bnd

end

function generate_bnds_ξ(ps, ed :: ETBDims)

    ps = copy(ps)

    ps[ (ps .* ( ed.μ  +  ed.ϵ)) .< 1] .= NaN

    Ψ_c  = ( ps .* (ed.μ + ed.ϵ) ).^0.5 .- 1   # Critical Ψ

    ξ_left_bnd  = ed.μ / ed.ν * ( (1 .- ed.ϵ .* Ψ_c ./ ed.μ) ./ (1 .+ abs.(Ψ_c)) - (1 .- Ψ_c / ed.μ) ./ ps )
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
    
    γ_left_bnd  = p_left_bnd  * ed.factors.p2γSv
    γ_right_bnd = p_right_bnd * ed.factors.p2γSv
    γs          = ps          * ed.factors.p2γSv


    ξ_left_bnd  .+= ed.ξ0
    ξ_right_bnd .+= ed.ξ0
    ξs          .+= ed.ξ0

    #println("ξ0 = ", ed.ξ0)
    #println("ξs + ξ0 = ", ξs .+ ed.ξ0)

    data[k] = Dict(
        "ξ_bnds" => (ξ_left_bnd, ξ_right_bnd),
        "p_bnds" => (p_left_bnd, p_right_bnd),
        "γ_bnds" => (γ_left_bnd, γ_right_bnd),
        "ps" => ps,
        "γs" => γs,
        "ξs" => ξs,
    )


end


regimes["etb_zatom"] = Dict(
        "label"     => "ETBM with \$100\\%\$ V",
        "label_pos" => (0.134, -0.4),
)

for k in keys(data)
    
    println(k)

    ξs = data[k]["ξs"]
    γs = data[k]["γs"]

    fixed_ξ = zeros(Float64, length(ξs), 3)
    fixed_γ = zeros(Float64, length(γs), 3)

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

plot_cases = ["etb_zatom", "standard"]
fcs = [ "none", "none", "none", "none"][end:-1:1]
ecs = [ "blue", "red", "green" ][end:-1:1]
hatches = [ "..", "//", "//", "//" ][end:-1:1]

fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel("\$\\gamma\$ [Sv]", fontsize=25)
ax.set_ylabel("\$\\xi\$", fontsize=25)
ax.grid()
ax.set_ylim([-1.6, 0.1])
ax.set_xlim([0.04, 0.2])

#ax.fill_betweenx(ξs, γ_left_bnd, γ_right_bnd, facecolor="blue", edgecolor="blue",       hatch="..", alpha=0.8, linewidth=1, zorder=10)#, label="Folding along fixed \$\\xi\$")
#ax.fill_between(ps * factor_p2γSv, ξ_left_bnd, ξ_right_bnd, facecolor="none",  edgecolor="orangered",  hatch="//", alpha=0.8, linewidth=1, zorder=10)#, label="Folding along fixed \$p\$")

#for (k, key) in enumerate(["standard", "etb_halfV", "etb_zatom", "etb_tenthV"])
for (k, key) in enumerate(plot_cases)

    println("Plotting regime = $key")

    regime = regimes[key]

    fixed_ξ = regime["fixed_ξ"]
    fixed_γ = regime["fixed_γ"]

    label = regime["label"]

    # creating polygons
    poly_ξ = createPoly("vertical",   fixed_ξ[:, 1], fixed_ξ[:, 2], fixed_ξ[:, 3])
    poly_γ = createPoly("horizontal", fixed_γ[:, 1], fixed_γ[:, 2], fixed_γ[:, 3])
    #merged_poly = shp_ops.unary_union([poly_ξ, poly_γ])
    merged_poly = shp_ops.unary_union([poly_γ, poly_ξ])
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

    #if length(ARGS) >= 4
    #    ax.set_title(ARGS[4])
    #end

    #ax.text(regime["label_pos"]..., regime["label"], size=12, ha="center", va="center", color=ecs[k])
end

#ax.scatter(P_adjusted..., s=20, color="black", zorder=99)
#ax.text(P_adjusted[1] + 0.01, P_adjusted[2], "\$P\$", va="center", ha="center", size=15)

ax.legend(loc="lower right")

fig.savefig(format("figures/regime_diagrams_comparison_xi.png"), dpi=300)

plt.show()
