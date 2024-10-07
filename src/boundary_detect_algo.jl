

function generate_bound_points(
    fixed_ξ :: AbstractArray{Float64, 2},
    fixed_γ :: AbstractArray{Float64, 2},
)
    
    println("Assume ξ and γ are sorted from small to large.")

    N_fixed_ξ = size(fixed_ξ, 1)
    N_fixed_γ = size(fixed_γ, 1)

    fulldata = Dict(
        :fixed_ξ => fixed_ξ,
        :fixed_γ => fixed_γ,
    )

    total_N = N_fixed_ξ + N_fixed_γ

    upper_bnds = zeros(Float64, total_N, 2)
    lower_bnds = zeros(Float64, total_N, 2)
    
    order = range(1, total_N) |> collect
    all_γ_top_branch = vcat(fixed_ξ[:, 2], fixed_γ[:, 1])
    all_γ_bot_branch = vcat(fixed_ξ[:, 3], fixed_γ[:, 1])
    order_within = vcat(collect(1:N_fixed_ξ), collect(1:N_fixed_γ))
    identities = vcat(
        repeat([:fixed_ξ], outer=N_fixed_ξ),
        repeat([:fixed_γ], outer=N_fixed_γ),
    )

    sort_idx_top_branch = sortperm(all_γ_top_branch)
    sort_idx_bot_branch = sortperm(all_γ_bot_branch)

    local pts = Dict(:top => [], :bot => [])

    for (branch, sort_idx) in Dict(
        :top => sort_idx_top_branch, 
        :bot => sort_idx_bot_branch,
    )
        

        for idx in sort_idx
            
            _identity = identities[idx]
            _order_within = order_within[idx]
            
            if _identity == :fixed_ξ
                if branch == :top
                    bnd_idx = 2
                else
                    bnd_idx = 3
                end
                _ξ = fulldata[_identity][_order_within, 1]
                _γ = fulldata[_identity][_order_within, bnd_idx]
            elseif _identity == :fixed_γ
                if branch == :top
                    bnd_idx = 3
                else
                    bnd_idx = 2
                end
                _γ = fulldata[_identity][_order_within, 1]
                _ξ = fulldata[_identity][_order_within, bnd_idx]
            end

            push!(pts[branch], [ _γ, _ξ, ])
            
        end
    end

    combined_pts = [ pts[:top]... , pts[:bot][end:-1:1]...]

    return combined_pts
end
