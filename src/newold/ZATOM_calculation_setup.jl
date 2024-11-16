
N = 1001

factor_PmE_mstommyear = 86400*365*1e3

γ0 = 1e6 # m^3/s


H = 4500.0
H_m = 60.0 # m
ξ = 0.5

ϕc = 40.0 |> deg2rad
ϕs = 10.0 |> deg2rad
ϕn = 70.0 |> deg2rad
δ = 20.0 |> deg2rad

Δλi = 50.0
Δλb =  5.0


S0 = 35.0 # PSU

Ts = 30.0
Tn = 5.0

a = 6.4e6

Δλe = Δλi + Δλb


function integrate(
    y_func,
    x0,
    x1,
    n ;
    if_cos :: Bool = true,
)

    xs = collect(range(x0, x1, length=n+1))
    dx = xs[2] - xs[1]
    y = y_func.(xs) .* ( (if_cos) ?  cos.(xs) : 1.0 )
    int_y = (y[1] + y[end] + 2 * sum(y[2:end-1])) * dx / 2

    return int_y
end

function my_tanh(ϕ)
    return tanh( (ϕ - ϕc) / δ ) 
end

Θϕc(ϕ) = ( ϕ > ϕc ) ? 1.0 : 0.0


balance_factor = - integrate(my_tanh, ϕs, ϕc, N) / integrate(my_tanh, ϕc, ϕn, N)
println("Balance factor: ", balance_factor)
σ(ϕ) = my_tanh(ϕ) * ( (ϕ > ϕc) ? balance_factor : 1.0 )

σ_pos(ϕ) = (σ(ϕ) + abs(σ(ϕ))) / 2.0
A0 = a^2 * (Δλe + Δλb) * integrate(σ_pos, ϕs, ϕn, N) # compute the total positive flux

function ew_factor(ϕ, ξ, side)
    
    if side == "w"
        factor = 1 - ξ * Θϕc(ϕ)
    elseif side == "e"
        factor = 1 + ξ * Δλb / Δλe * Θϕc(ϕ)
    else
        throw(ErrorException("Unknown side: $side"))
    end
    return factor
end

PmE_w(ϕ) = γ0 * ew_factor(ϕ, ξ, "w") * σ(ϕ) / A0
PmE_e(ϕ) = γ0 * ew_factor(ϕ, ξ, "e") * σ(ϕ) / A0
PmE_ξ0(ϕ) = γ0 * ew_factor(ϕ, 0, "w") * σ(ϕ) / A0
FWF(ϕ) = ( PmE_w(ϕ) * Δλb + PmE_e(ϕ) * Δλe ) * a * cos(ϕ)
T_sfc(ϕ) = Tn + (Ts - Tn) * 0.5 * (1 + cos( (ϕ - ϕs) / (ϕn - ϕs) * π) ) 

println("Verify the total flux (Sv): ", a * integrate(FWF, ϕs, ϕn, N; if_cos=false) / 1e6)
println("Verify the total pos flux (Sv): ", a * integrate(FWF, ϕc, ϕn, N; if_cos=false) / 1e6)



println("Compute various integration constants.")

ϕ_1 = 40 |> deg2rad
ϕ_2 = 68 |> deg2rad
k = 1
tilde_k_Θ = 2 / (k * π) * ( 1 - cos(H_m / H * k * π) )

σΘϕc(ϕ) = σ(ϕ) * Θϕc(ϕ)

flat(ϕ) = 1.0
mean_σΘϕc = integrate(σΘϕc, ϕ_1, ϕ_2, N) / integrate(flat, ϕ_1, ϕ_2, N)

println("k = ", k)
println("tilde_k_Θ = ", tilde_k_Θ)
println("mean_σΘϕc = ", mean_σΘϕc)
println("A0 = ", A0)
println("H_m = ", H_m)

α_S = 0.7 # kg/m^3/PSU
g = 10.0 # m/s^2
ρ0 = 1e3 # kg/m^3

println("tilde_k_Θ * mean_σΘϕc * S0 * (1+Δλb/Δλe) * α_S * g / ρ0 / (H_m * A0)= ", tilde_k_Θ * mean_σΘϕc * S0 * (1+Δλb/Δλe) * α_S * g / ρ0 / (A0 * H_m) )






