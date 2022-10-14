using Formatting
using Roots
using PyCall

shp_geo = pyimport("shapely.geometry")
shp_ops = pyimport("shapely.ops")
plt_patches = pyimport("matplotlib.patches")

include("ZATOM_regimes.jl")

regimes["standard"]["label_pos"] = (0.11, -1.1)
regimes["standard"]["label"] = "ZATOM"

ξ_bnd = -0.9
Ae = 3.0
Aw = 1.0
AH = Ae + Aw
AL = AH * 1
a = 0.475 # ( 1 - ξ_bnd ) / (1 + Ae / Aw)
γ_hyds = collect(range(0, 0.5e6, length=1001))
γ_locs = [0.001, 0.005, 0.01, 0.015] * 1e6 #collect(range(0, 0.01e6, length=11))[2:end]

γ_hyds_marks = collect(range(0, 0.5e6, length=11))

function tranform_coor(γ_hyd, γ_loc)
    r = γ_loc / γ_hyd
    return (
        γ_hyd + AL / (Ae + AL) * γ_loc,
        1 - ( 1 + Ae/Aw ) * (γ_loc + a * γ_hyd) / (γ_hyd + γ_loc * (AL / (Ae + AL))),
    )
    
end



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

plot_cases = ["standard", ]
fcs =        [ "none", "none", "none", "none", "none", "none"]
ecs =        [ "black", "red", "blue", "green", "orange", "gray"]
hatches =    [ "..", "..", "\\\\", "//", "||", "||"]
zorders =    [ 10, 1, 5, 4, 3, 2 ]
fig, axes = plt.subplots(1, 2, constrained_layout=true, figsize=(10,4))


for l = 1:2

    ax = axes[l]
    ax.set_xlabel("\$\\gamma\$ [Sv]", fontsize=25)
    ax.set_ylabel("\$\\xi\$", fontsize=25)
    ax.grid(alpha=0.5)

    for (k, key) in enumerate(plot_cases)

        println("Plotting regime = $key")

        regime = regimes[key]

        fixed_ξ = regime["fixed_ξ"]
        fixed_γ = regime["fixed_γ"]

        label = regime["label"]

        # creating polygons
        poly_ξ = createPoly("vertical",   fixed_ξ[:, 1], fixed_ξ[:, 2], fixed_ξ[:, 3])
        poly_γ = createPoly("horizontal", fixed_γ[:, 1], fixed_γ[:, 2], fixed_γ[:, 3])


        merged_poly = shp_ops.unary_union([poly_γ, poly_ξ])
        merged_poly = shpPoly2MatplotPoly(merged_poly, Dict(
            :ec     => ecs[k],
            :fc     => fcs[k],
            :alpha  => 1.0,
            :zorder => zorders[k],
            :hatch  => hatches[k],
            :label  => regime["label"],
        ))
        
        ax.add_patch(merged_poly)

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

    for (i, γ_loc) in enumerate(γ_locs)

        xs = []
        ys = []

        xys = []
        for (k, γ_hyd) in enumerate(γ_hyds)
            x, y = tranform_coor(γ_hyd, γ_loc)
            push!(xs, x)
            push!(ys, y)
        end
        xs ./= 1e6
        ax.plot(xs, ys, label=format("\$\\gamma' = {:.3f} \\, \\mathrm{{Sv}}\$", γ_loc/1e6))


        xs = []
        ys = []
        xys = []
        for (k, γ_hyd) in enumerate(γ_hyds_marks)
            x, y = tranform_coor(γ_hyd, γ_loc)
            push!(xs, x)
            push!(ys, y)
        end
        xs ./= 1e6
        ax.scatter(xs, ys)


        #println(xs)
        #println(ys)
    end
 

    if l == 1

        ax.set_title("(a)")

        ax.legend(loc="lower right")
        ax.set_ylim([-5, -0.6])
        ax.set_xlim([-0.01, 0.15])

        ax.text(0.03, -4.5, format("\$\\overline{{\\gamma}} = {:d} \\, \\mathrm{{Sv}}\$", γ_hyds_marks[1] / 1e6), ha="center", va="center")
        ax.text(0.0655, -1.86, format("\$\\overline{{\\gamma}} = {:.2f} \\, \\mathrm{{Sv}}\$", γ_hyds_marks[2] / 1e6), ha="center", va="center")
        ax.text(0.1157, -1.55, format("\$\\overline{{\\gamma}} = {:.1f} \\, \\mathrm{{Sv}}\$", γ_hyds_marks[3] / 1e6), ha="center", va="center")
    elseif l==2
        ax.set_title("(b)")
        ax.set_ylim([-1.6, -0.7])
        ax.set_xlim([0.04, 0.12])
        ax.text(0.0748, -1.551, format("\$\\overline{{\\gamma}} = {:.2f} \\, \\mathrm{{Sv}}\$", γ_hyds_marks[2] / 1e6), ha="center", va="center")
        ax.text(0.1093, -1.371, format("\$\\overline{{\\gamma}} = {:.1f} \\, \\mathrm{{Sv}}\$", γ_hyds_marks[3] / 1e6), ha="center", va="center")

    end

end 


fig.savefig(format("figures/figure-approx_dijkstra.png"), dpi=300)

plt.show()
