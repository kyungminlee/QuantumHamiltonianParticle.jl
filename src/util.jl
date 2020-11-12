export tuplelength

tuplelength(::Type{<:NTuple{N, Any}}) where N = N
tuplelength(::NTuple{N, Any}) where N = N

tupleadd(l::T, r::T) where {T<:Tuple} = l .+ r
tuplesubtract(l::T, r::T) where {T<:Tuple} = l .- r

tuplezero(::Type{T}) where {T<:Tuple} = ((zero(S) for S in T.parameters)...,)
tupleone(::Type{T}) where {T<:Tuple} = ((one(S) for S in T.parameters)...,)

tuplezero(::T) where {T<:Tuple} = ((zero(S) for S in T.parameters)...,)
tupleone(::T) where {T<:Tuple} = ((one(S) for S in T.parameters)...,)

export isparityodd

# in utility
function isparityodd(arr::AbstractVector{T})::Bool where T
    n = length(arr)
    parity = false
    for i in n-1:-1:1, j in i+1:n
        if arr[j] < arr[i]
            parity = !parity
        end
    end
    return parity
end
