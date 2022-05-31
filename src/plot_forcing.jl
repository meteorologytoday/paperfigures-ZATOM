using Formatting

N = 1001

ϕc = 40.0 |> deg2rad
ϕs = 10.0 |> deg2rad
ϕn = 70.0 |> deg2rad
δ = 20.0 |> deg2rad

Ts = 30.0
Tn = 5.0


function integrate(y_func, x0, x1, n)
    xs = collect(range(x0, x1, length=n+1))
    dx = xs[2] - xs[1]
    y = y_func.(xs) .* cos.(xs)
    int_y = (y[1] + y[end] + 2 * sum(y[2:end-1])) * dx

    return int_y
end

function my_tanh(ϕ)
    return tanh( (ϕ - ϕc) / δ ) 
end

balance_factor = - integrate(my_tanh, ϕs, ϕc, N) / integrate(my_tanh, ϕc, ϕn, N)
println("Balance factor: ", balance_factor)
σ_unorm(ϕ) = my_tanh(ϕ) * ( (ϕ > ϕc) ? balance_factor : 1.0 )
σ_unorm_pos(ϕ) = max(σ_unorm(ϕ), 0.0)
σ0 = integrate(σ_unorm_pos, ϕs, ϕn, N)

σ(ϕ) = σ_unorm(ϕ) / σ0
σ_pos(ϕ) = max(σ(ϕ), 0.0)

T_sfc(ϕ) = Tn + (Ts - Tn) * 0.5 * (1 + cos( (ϕ - ϕs) / (ϕn - ϕs) * π) ) 

println("Verify the total flux: ", integrate(σ, ϕs, ϕn, N))
println("Verify the total pos flux: ", integrate(σ_pos, ϕs, ϕn, N))

println("Loading PyPlot")
using PyPlot
plt = PyPlot
plt.ion()
println("Done")

fig, ax = plt.subplots(1, 1, figsize=(6,4), constrained_layout=true)

ax_twinx = ax.twinx()

ϕ = collect(range(ϕs, ϕn, length=N))

ln_T_sfc = ax.plot(rad2deg.(ϕ), T_sfc.(ϕ), "k-")
ln_σ = ax_twinx.plot(rad2deg.(ϕ), σ.(ϕ), "r-")

ax.set_xlabel("\$\\phi\$ [ \${}^{\\circ}\$N ]")

ax.set_ylabel("\$ T_{\\mathrm{sfc}} \$ [ \${}^{\\circ}\\mathrm{C} \$ ]")
ax_twinx.set_ylabel("\$\\sigma \$ [ \$\\mathrm{m} / \\mathrm{s} \$ ]", color="red")
ax_twinx.spines["right"].set_color("red")


ax_twinx.set_yticks([5, 10, 15, 20, 25, 30])
ax_twinx.set_yticks([-1, -0.5, 0, 0.5, 1.0, 1.5])
ax.set_ylim([2.5, 32.5])
ax_twinx.set_ylim([-1.25, 1.75])
ax_twinx.yaxis.label.set_color("red")
ax_twinx.tick_params(axis="y", colors="red")
ax.grid()


fig.savefig("figure-forcing.png", dpi=300)

sleep(500)
