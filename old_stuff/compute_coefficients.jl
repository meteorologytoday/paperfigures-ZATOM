using Formatting

# The following two boundaries are used to
# compute the averaged values for diagnostic
# equation
avg_ϕs = 40.0 |> deg2rad
avg_ϕn = 70.0 |> deg2rad
H_S = 60.0
α_S = 0.7  # kg/m^3/PSU
S0  = 35.0 # PSU
ρ0  = 1000.0 # kg/m^3
Re = 6.4e6
Ω = 7.3e-5
g0 = 10.0  # m/s^2
f_co(ϕ) = 2 * Ω * sin(ϕ)

N = 1001

ϕc = 40.0 |> deg2rad
ϕs = 10.0 |> deg2rad
ϕn = 70.0 |> deg2rad
δ = 20.0 |> deg2rad

Δλi = 10.0
Δλb =  5.0
D = Δλi / (Δλi + Δλb)
Λ = Δλb / (Δλi + Δλb)
ξ = 0.2


Ts = 30.0
Tn = 5.0

H = 4500.0

ψ1(z) = - sin(π*z/H)
function avg_z(z_func, n)
    
    zs = collect(range(-H, 0, length=n+1))
    dz = zs[2] - zs[1]
    z = z_func.(zs) .* ψ1.(zs)
    result = (z[1] + z[end] + 2 * sum(z[2:end-1])) * dz / 2
    result *= 2 / H

    return result
   

end

function integrate(y_func, x0, x1, n; avg :: Bool = false)
    xs = collect(range(x0, x1, length=n+1))
    dx = xs[2] - xs[1]
    y = y_func.(xs) .* cos.(xs)
    int_y = (y[1] + y[end] + 2 * sum(y[2:end-1])) * dx / 2

    if avg
        int_y /= (sin(x1) - sin(x0))
    end

    return int_y
end

function my_tanh(ϕ)
    return tanh( (ϕ - ϕc) / δ ) 
end

balance_factor = - integrate(my_tanh, ϕs, ϕc, N) / integrate(my_tanh, ϕc, ϕn, N)
println("Balance factor: ", balance_factor)
σ_unorm(ϕ) = my_tanh(ϕ) * ( (ϕ > ϕc) ? balance_factor : 1.0 )
σ_unorm_pos(ϕ) = max(σ_unorm(ϕ), 0.0)

σ0 = integrate(σ_unorm_pos, ϕs, ϕn, N) * Re^2 * deg2rad( 2 * Δλb + Δλb )

σ(ϕ) = σ_unorm(ϕ) / σ0
σ_pos(ϕ) = max(σ(ϕ), 0.0)

function σ_ξ(ϕ, ξ; side)

    local _σ

    _σ = σ(ϕ)
   
    if ϕ >= ϕc

        if side == "w"
            _σ *= 1 - ξ
        elseif side == "e"
            _σ *= 1 + ξ * D
        else
            throw(ErrorException("Unknown side: $side"))
        end
    end

    return _σ

end

σ_w(ϕ) = σ_ξ(ϕ, ξ; side="w")
σ_e(ϕ) = σ_ξ(ϕ, ξ; side="e")


T_sfc(ϕ) = Tn + (Ts - Tn) * 0.5 * (1 + cos( (ϕ - ϕs) / (ϕn - ϕs) * π) ) 

println("Verify the total flux: ", integrate(σ, ϕs, ϕn, N))
println("Verify the total pos flux: ", integrate(σ_pos, ϕs, ϕn, N))

# Test
println("Test: (2/H) ∫ sin^2(πz/H) dz = ", avg_z(ψ1, 1000), ". Expected value = 1")

# Compute <f>
f_avg = 4 * Ω / π * (sin(avg_ϕn) - sin(avg_ϕs))
println("<f> = ", f_avg)
println("A = (π / H)^2 * <f> = ", (π / H)^2 * f_avg)

# Compute <σ Θ_ϕc Θ_HS>

## First: compute avg_y of σ Θ_ϕc
Θ_ϕc(ϕ) = (ϕ > ϕc) ? 1.0 : 0
avg_σ_Θ_ϕc = integrate((ϕ,)-> σ(ϕ) .* Θ_ϕc(ϕ), avg_ϕs, avg_ϕn, 1000; avg = true)
println("<σ Θ_ϕc> = ", avg_σ_Θ_ϕc)

## Second: compute avg_z of Θ_HS
Θ_HS(z) = (z > - H_S) ? 1.0 : 0
avg_Θ_HS = avg_z(Θ_HS, 10000)
println("<Θ_HS> = ", avg_Θ_HS)

## Final: compute <σ Θ_ϕc Θ_HS>
avg_σ_Θ_ϕc_Θ_HS = avg_σ_Θ_ϕc * avg_Θ_HS
println("<σ Θ_ϕc Θ_HS> = ", avg_σ_Θ_ϕc_Θ_HS)


println("Therefore: ")
println("α_S g / ρ0 (1 + Λ) S0 / H_S = ", α_S * g0 / ρ0 * (1 + Λ) * S0 / H_S)
println("α_S g / ρ0 (1 + Λ) S0 <σ Θ_ϕc Θ_HS> / H_S = ", avg_σ_Θ_ϕc_Θ_HS * α_S * g0 / ρ0 * (1 + Λ) * S0 / H_S)






