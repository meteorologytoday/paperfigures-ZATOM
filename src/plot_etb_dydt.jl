using Formatting
println("Loading ArgParse and PyPlot...")
using ArgParse
using PyPlot
plt = PyPlot
println("Done")

using Roots

s = ArgParseSettings()
@add_arg_table! s begin

    "--title"
        help = "Title of the figure"
        arg_type = String
        default = ""

    "--eps"
        help = "Value of epsilon ϵ"
        arg_type = Float64
        required = true

    "--find-zeros-init-guess"
        help = "Value of initial guess of on, middle, and off state"
        nargs = 3
        arg_type = Float64
        default = [0.0, 0.5, 1.1]

    "--output"
        help = "Output file name."
        arg_type = String
        default = ""

       
end

parsed = parse_args(ARGS, s)

# ===== Figure varying ξ =====

p = 1.2
ξ = 0.0
μ = 5.0
ν = 1.0
ϵ = parsed["eps"]

y_max = (ϵ == 0) ? 1.5 : 1.0 / ϵ
sampling_density = 500 # pts per unit length


pos(x) = (x >= 0) ? x : 0.0
Ψ(y, p, ξ, μ, ν, ϵ)         = ( μ * (1 - y) + ν * p * ξ ) / ( 1 - pos(ϵ * y) )
_dydt(y, p, ξ, μ, ν, ϵ)     = p - ( 1 + abs( Ψ(y, p, ξ, μ, ν, ϵ) ) ) * y
_y_kink(p, ξ, μ, ν, ϵ) = 1 + ν * p * ξ / μ 

y_interval = [-0.1, y_max]
ys = collect(range(y_interval[1], y_interval[2], length=round(Integer, (y_interval[2] - y_interval[1]) * sampling_density)))
dydt(p, ξ, μ, ν, ϵ) = [ _dydt(y, p, ξ, μ, ν, ϵ) for y in ys ]

example_coe = (p, ξ, μ, ν, ϵ)


fig, ax = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=true)
for ξ in [-1.0, -2.0, 1.0, 0.0]

    if ξ == 0.0
        ax.plot(ys, dydt(p, ξ, μ, ν, ϵ), "-", color="#000000", zorder=10)
    else
        ax.plot(ys, dydt(p, ξ, μ, ν, ϵ), "-", color="#cccccc", zorder=1)
    end
    y_kink = _y_kink(p, ξ, μ, ν, ϵ)
    ax.text(y_kink + ((ξ==0) ? 0.1 : 0), _dydt(y_kink, p, ξ, μ, ν, ϵ) + 0.1, "\$ \\tilde{\\xi} = $(format("{:d}", ξ)) \$", color= ( (ξ == 0) ? "black" : "#aaaaaa" ), ha="center", va="bottom")
end

ax.plot([-5, 5], [0, 0], color="black", linewidth=1)
ax.plot([0, 0], [-5, 5], color="black", linewidth=1)

x0 = 0
x1 = 1 + ν * p * ξ / μ
y0 = _dydt(x0, example_coe...)
y1 = _dydt(x1, example_coe...)
ax.scatter([x0, x1], [y0, y1], s=50, marker="o", color="red", zorder=20)
#ax.annotate("\$ A_1 = \\left(0, p\\right) \$", xy=(x0, y0),  xycoords="data",
ax.annotate("\$ A_1 \$", xy=(x0, y0),  xycoords="data",
            xytext=(0.3, 1.25), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=15,
)

#ax.annotate("\$ A_2 = \\left(1 + \\frac{\\nu p \\xi}{\\mu}, \\left(1 - \\frac{\\nu}{\\mu} \\xi \\right) p - 1\\right) \$", xy=(x1, y1),  xycoords="data",
ax.annotate("\$ A_2 \$", xy=(x1, y1),  xycoords="data",
            xytext=(1, 1), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=15,
)

#ax.text(x0 - 0.05, y0, "\$ \\left(0, p\\right) \$", va="center", ha="center", fontsize=15)
#ax.text(x1, y1 + 0.25, "\$ \\left(1 + \\frac{\\nu}{\\mu} \\xi \\right) p - 1 \$", va="center", ha="center", fontsize=15)

slope = ( _dydt(0, example_coe...) - _dydt(1, example_coe...)) / (0 - 1)
ax.plot([0, 2], [_dydt(0, example_coe...), _dydt(0, example_coe...) + 2 * slope], color="red", linestyle="dashed")

# Dot steady states

f(y) = _dydt(y, example_coe...)

y_ss    = [ find_zero(f,  init_pt) for init_pt in parsed["find-zeros-init-guess"] ]
dydt_ss = f.(y_ss)
ax.scatter(y_ss, dydt_ss, s=50, marker="x", color="red", zorder=20)

ax.plot([y_ss[1], y_ss[1]], [ 0, -0.5 ], color="#888888", linestyle="dashed")
ax.plot([y_ss[3], y_ss[3]], [ 0, -0.5 ], color="#888888", linestyle="dashed")
ax.text(y_ss[1], -0.5, "on", color="red", ha="center", va="top")
ax.text(y_ss[3], -0.5, "off", color="red", ha="center", va="top")


# Dot the intersected point
#ax.scatter([p,], [0,], s=10, marker="o", color="black", zorder=50)
#ax.text(p, - 0.2, "\$\\left(p, 0\\right)\$", color="black", ha="center", va="top", fontsize=11)

# Dot the lowest point on the parabola

interval_of_interest = [0, x1]
ys    = collect(range(interval_of_interest..., length=500))
dydts = [ _dydt(y, example_coe...) for y in ys ]
low_idx = argmin(dydts)

y_low = ys[low_idx]
dydt_low = dydts[low_idx]
ax.scatter([y_low], [dydt_low], s=50, marker="o", color="red", zorder=20)

#ax.annotate("\$ A_3 = \\left( \\frac{1 + \\mu + \\nu p \\xi}{2 \\mu} , p - \\frac{\\left(1 + \\mu + \\nu p \\xi\\right)^2}{4 \\mu}\\right) \$", xy=(y_low, dydt_low),  xycoords="data",
ax.annotate("\$ A_3 \$", xy=(y_low, dydt_low),  xycoords="data",
            xytext=(y_low, dydt_low - 0.7), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=15
)


# Dot the intersected point
ax.scatter(p, 0, s=50, marker="o", color="red", zorder=20)
ax.annotate("\$ A_4 \$", xy=(p, 0),  xycoords="data",
            xytext=(p, -0.7), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=15
)


ax.set_title(parsed["title"])


ax.set_ylabel("\$\\mathrm{d} y / \\mathrm{d}\\tau\$")
ax.set_xlabel("\$y\$")

ax.grid()
ax.set_ylim([-2, 1.5])
ax.set_xlim([-0.1, 1.4])

if parsed["output"] != ""
    fig.savefig(parsed["output"], dpi=300)
end

#plt.show()


