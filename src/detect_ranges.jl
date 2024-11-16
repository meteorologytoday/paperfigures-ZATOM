





function detectRanges(arr, include_prev=false)

    vals = []
    rngs = []

    N = length(arr)

    new_rng = true
    beg_i = 1
    detect_val = 0.0

    extra_include = 0
    if include_prev
        extra_include = -1
    end

    for i in 1:N
    
        if new_rng
            new_rng = false
            detect_val = arr[beg_i]
        else
            if arr[i] == detect_val
                # do nothing
            else
                new_rng = true
            end
        end

        if new_rng || i == N

            push!(vals, detect_val)

            if new_rng
                push!(rngs, beg_i+extra_include:i-1)
            elseif i == N
                push!(rngs, beg_i+extra_include:i)
            end

            beg_i = i

        end

    end

    return vals, rngs
end

