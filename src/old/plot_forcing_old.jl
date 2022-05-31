using Formatting

ϕc = 40.0 |> deg2rad
ϕs = 10.0 |> deg2rad
ϕn = 70.0 |> deg2rad
ξ = 0.0
D = 5.0 / 20.0

function integrate(y_func, x0, x1, n)
    xs = collect(range(x0, x1, length=n+1))
    dx = xs[2] - xs[1]
    y = y_func.(xs)
    int_y = (y[1] + y[end] + 2 * y[2:end-1]) * dx

    return int_y
end

function my_tanh(ϕ)
    return tanh( (ϕ - ϕc) / δ ) 
end

G_B_factor = - integrate(my_tanh, ϕs, ϕc) / integrate(my_tanh, ϕc, ϕn)

G_A(ϕ) = 1.0
G_B(ϕ) = ( ϕ > ϕc ) ? 1.0 : 1.0


function μ(ϕ, form)
    if form == :A
        G = G_A
    elseif form == :B
        G = G_B
    end
    return my_tanh(ϕ) * G(ϕ)
end

function r(ϕ, Ω)
    if ϕ > ϕc && Ω == :w
        return 1.0 + ξ
    elseif ϕ > ϕc && Ω == :e
        return 1.0 - ξ * D
    else
        return 1.0
    end
end

println("Loading PyPlot")
using PyPlot
plt = PyPlot
plt.ion()
println("Done")

fig, ax = plt.subplots(2, 1, constrained_layout=true)
N = 101
ϕ = collect(range(ϕs, ϕn, length=101))
ax.plot(ϕ, r(ϕ))



ax.set_xlabel("\$\\phi\$")
ax.set_ylabel("\$\\Theta\$")
ax.grid()

ax.legend()
fig.savefig("figure-forcings.png", dpi=300)

sleep(500)
