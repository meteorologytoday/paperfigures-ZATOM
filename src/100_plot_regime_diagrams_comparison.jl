using Formatting
using Roots
using PyCall

shp_geo = pyimport("shapely.geometry")
shp_ops = pyimport("shapely.ops")
plt_patches = pyimport("matplotlib.patches")

include("ZATOM_regimes.jl")
include("boundary_detect_algo.jl")

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

    function ETBDims(
        t_d,
        t_r,
        V,
        δT_star,
        γ_cusp,
        ξ_cusp,
        Δξ_fold,
        Δξ_arc,
        ξ0;
        ϵ_zero :: Bool = false,
    )

        factor_γ2p = (2 * α_S * t_d * S0) / (α_T * δT_star * V)
        factor_p2γ = factor_γ2p^(-1)
        factor_p2γSv = factor_p2γ / 1e6
        
        p_cusp = γ_cusp * factor_γ2p

        # Method 1
        #=
        r = Δξ_fold / Δξ_arc
        μ = ( ( 1 / p_cusp - 1) / r)^0.5
        ν = μ^2 / Δξ_arc
        ϵ = Δξ_fold * ν + 1 - μ
        =#


        # Method 2

        if ϵ_zero == false
            μ = Δξ_arc / (ξ_cusp - ξ0) * ( 1 - p_cusp^(-1))
            ν = μ^2 / Δξ_arc
            ϵ = 1 / p_cusp - μ
        else
            μ = p_cusp^(-1)
            ν = μ * (1 - μ) / (ξ_cusp - ξ0)
            ϵ = 0.0 
        end

        println("p_cusp = $(p_cusp)")
        println("V = $(V/1e15) * 1e15 m^3")
        println("t_d = $(t_d / 86400/365) yr.")
        println("factor_γ2p = $factor_γ2p")
        println("μ = $μ")
        println("ν = $ν")
        println("ϵ = $ϵ")
        
        println("Upper limit for ξ = μ/ν = $(μ / ν)")

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


ξs_all = collect(range( -5,  0.5,  length=101))
ps_all = collect(range( -0.1, 10.0,  length=1001))

ϕn = 70.0
ϕs = 10.0
H = 4500.0
a = 6.4e6
α_T = 2e-3
α_S = 7e-3

# divided by 2 is important
full_V   = ((60 / 360) * 2π * H * a^2 * (sin(deg2rad(ϕn)) - sin(deg2rad(ϕs))) / 2)
V        = full_V  *  (500/H) * (5/60)

println("full_V = ", full_V)
println("V / full_V = ", V / full_V)

δT_star = 25.0
S0 = 35.0

L = a * deg2rad(ϕn - ϕs)
Kh = 4e4

t_d = L^2 / Kh 
t_a = V / t_d 
t_r = 10 * 86400.0

Q = t_d / t_r

ξ0 = -4.3
cusp_pt = (γ=0.025e6, ξ=-5.7498)
ξ_u = -3.7

γ_cusp  = cusp_pt.γ
ξ_cusp  = cusp_pt.ξ


Δξ_fold = 0.0
Δξ_arc  = ξ_u - ξ_cusp

println("Δξ_arc = ", Δξ_arc)



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
        ξ_cusp,
        Δξ_fold,
        Δξ_arc,
        ξ0;
    ),

    "etb_zatom_eps0" => ETBDims(
        t_d,
        t_r,
        V,
        δT_star,
        γ_cusp,
        ξ_cusp,
        Δξ_fold,
        Δξ_arc,
        ξ0;
        ϵ_zero = true,
    ),

    "etb_zatom_V0.8" => ETBDims(
        t_d,
        t_r,
        0.8 * V,
        δT_star,
        γ_cusp,
        ξ_cusp,
        Δξ_fold,
        Δξ_arc,
        ξ0;
    ),

    "etb_zatom_V1.2" => ETBDims(
        t_d,
        t_r,
        1.2 * V,
        δT_star,
        γ_cusp,
        ξ_cusp,
        Δξ_fold,
        Δξ_arc,
        ξ0;
    ),

    "etb_zatom_xi" => ETBDims(
        t_d,
        t_r,
        V,
        δT_star,
        γ_cusp,
        ξ_cusp,
        Δξ_fold,
        Δξ_arc,
        -4;
    ),

    "etb_zatom_xi_0" => ETBDims(
        t_d,
        t_r,
        V,
        δT_star,
        γ_cusp,
        ξ_cusp,
        Δξ_fold,
        Δξ_arc,
        ξ0 + 0.5;
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

    println("## Computing boundaries for case: ", k)

    ξs, p_left_bnd, p_right_bnd = generate_bnds_p(ξs_all, ed)
    ps, ξ_left_bnd, ξ_right_bnd = generate_bnds_ξ(ps_all, ed)


    #println("ps_all = ", ps_all)
    #println("ξ_left_bnd = ", ξ_left_bnd)
    #println("ξ_right_bnd = ", ξ_right_bnd)

    γ_left_bnd  = p_left_bnd  * ed.factors.p2γSv
    γ_right_bnd = p_right_bnd * ed.factors.p2γSv
    γs          = ps          * ed.factors.p2γSv


    ξ_left_bnd  .+= ed.ξ0
    ξ_right_bnd .+= ed.ξ0
    ξs          .+= ed.ξ0

    #println("γs = ", γs)
    #println("ξ_left_bnd = ", ξ_left_bnd)
    #println("ξ_right_bnd = ", ξ_right_bnd)
        
    println("Upper limit for ξ = μ/ν + ξ0 = $(ed.μ / ed.ν + ed.ξ0)")


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
        "label"     => "ETBM",
)

regimes["etb_zatom_eps0"] = Dict(
        "label"     => "ETBM with \$\\epsilon = 0\$",
)

regimes["etb_zatom_V0.8"] = Dict(
        "label"     => "ETBM with \$0.8 \\, V\$",
)

regimes["etb_zatom_V1.2"] = Dict(
        "label"     => "ETBM with \$1.2 \\, V\$",
)

regimes["etb_zatom_xi"] = Dict(
        "label"     => "ETBM with \$ \\xi_0 = -4.0 \$",
)

regimes["etb_zatom_xi_0"] = Dict(
        "label"     => "ETBM with \$ \\xi_0 \$ increase by \$0.5\$ ",
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

plot_cases = ["standard", "etb_zatom", "etb_zatom_eps0", "etb_zatom_V0.8", "etb_zatom_xi_0"]
fcs =        [ "none", "none", "none", "none", "none", "none"]
ecs =        [ "black", "red", "blue", "green", "orange", "gray"]
hatches =    [ "..", "..", "\\\\", "//", "||", "||"]
zorders =    [ 10, 1, 5, 4, 3, 2 ]
fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel("\$\\gamma\$ [Sv]", fontsize=25)
ax.set_ylabel("\$\\xi\$", fontsize=25)
ax.grid(alpha=0.5)

ax.set_ylim([-6, -3])
ax.set_xlim([0.0, 0.15])

for (k, key) in enumerate(plot_cases)

    println("Plotting regime = $key")

    regime = regimes[key]

    fixed_ξ = regime["fixed_ξ"]
    fixed_γ = regime["fixed_γ"]

    label = regime["label"]


    if key == "standard"

        pts = generate_bound_points(fixed_ξ, fixed_γ)

        println(pts)

        merged_poly = shp_geo.Polygon(pts)
        
        
    else
        # creating polygons
        poly_ξ = createPoly("vertical",   fixed_ξ[:, 1], fixed_ξ[:, 2], fixed_ξ[:, 3])
        poly_γ = createPoly("horizontal", fixed_γ[:, 1], fixed_γ[:, 2], fixed_γ[:, 3])


        merged_poly = shp_ops.unary_union([poly_γ, poly_ξ])

    end

    merged_poly = shpPoly2MatplotPoly(merged_poly, Dict(
        :ec     => ecs[k],
        :fc     => fcs[k],
        :alpha  => 1.0,
        :zorder => zorders[k],
        :hatch  => hatches[k],
        :label  => regime["label"],
    ))
 
    ax.add_patch(merged_poly)
    println("Merged_poly: ", merged_poly)

    # Plot the lines used to determine the boundaries
    if key == "standard"
    
        for i in 1:size(fixed_γ, 1)
            _γ       = fixed_γ[i, 1]
            _ξ_lower = fixed_γ[i, 2]
            _ξ_upper = fixed_γ[i, 3]
            ax.plot(
                [_γ, _γ, ],
                [_ξ_lower, _ξ_upper, ],
                linestyle="solid",
                color=ecs[k],
                linewidth=2,
                marker="o",
                markersize=5,
            )
        end

        for i in 1:size(fixed_ξ, 1)
            _ξ       = fixed_ξ[i, 1]
            _γ_lower = fixed_ξ[i, 2]
            _γ_upper = fixed_ξ[i, 3]
            ax.plot(
                [_γ_lower, _γ_upper, ],
                [_ξ, _ξ, ],
                linestyle="solid",
                color=ecs[k],
                linewidth=2,
                marker="o",
                markersize=5,

            )
        end 

    end

    # Plot the little rectangle to indicate the regime
    # is extending.
    if key == "standard"
        
        start_γ = 0.14119
        _ξ_bnd = [-0.9500, -0.6814]
        wedge_spacing_γ = 0.0015
        wedge_width_γ   = 0.0015
        slope = 1
        cnt = 3
        
        for i = 1:cnt
            dγ = (wedge_spacing_γ + wedge_width_γ) * (i-1)
            _start_γ = start_γ + dγ
            _ξ_bnd_shifted = _ξ_bnd .+ slope * dγ
            poly = shp_geo.Polygon([
                [_start_γ, _ξ_bnd_shifted[1]],
                [_start_γ + wedge_width_γ, _ξ_bnd_shifted[1] + wedge_width_γ * slope * 0],
                [_start_γ + wedge_width_γ, _ξ_bnd_shifted[2] + wedge_width_γ * slope * 0],
                [_start_γ, _ξ_bnd_shifted[2]]
            ])
            
            poly = shpPoly2MatplotPoly(poly, Dict(
                :ec     => ecs[k],
                :fc     => fcs[k],
                :alpha  => 1.0,
                :zorder => zorders[k],
                :hatch  => hatches[k],
            ))
            
            ax.add_patch(poly)

        end


    end

end

ax.legend(loc="lower right")

fig.savefig(format("figures/regime_diagrams_comparison.svg"), dpi=300)

plt.show()
