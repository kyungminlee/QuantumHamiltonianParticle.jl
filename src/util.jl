export tuplelength



tuplelength(::Type{<:NTuple{N, Any}}) where N = N
tuplelength(::NTuple{N, Any}) where N = N

tuplezero(::Type{<:NTuple{N, Any}}) where N = ([0 for i in 1:N]...,)
tuplezero(::NTuple{N, Any}) where N = tuplezero(NTuple{N, Any})


tupleadd(x::NTuple{N, Any}, y::NTuple{N, Any}) where {N} = x .+ y

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
