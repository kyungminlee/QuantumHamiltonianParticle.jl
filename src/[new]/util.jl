export tuplelength

tuplelength(::Type{<:NTuple{N, Any}}) where N = N
tuplelength(::NTuple{N, Any}) where N = N

export isparityodd

# in utility
function isparityodd(arr::AbstractVector{T})::Bool where T
    n = length(arr)
    parity = true
    for i in n-1:-1:1, j in i+1:n
        if arr[j] < arr[i]
            parity = !parity
        end
    end
    return parity
end
