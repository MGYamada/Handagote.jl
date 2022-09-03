struct DualArray{T <: ReComp, N} <: DenseArray{Dual{T}, N}
    value::Array{T, N}
    epsilon::Array{T, N}
    function DualArray{T, N}(A::AbstractArray{T, N}, B::AbstractArray{T, N}) where {T <: ReComp, N}
        @assert size(A) == size(B)
        new{T, N}(convert(Array, A), convert(Array, B))
    end
end

DualArray{T, N}(A::AbstractArray{T, N}) where {T <: ReComp, N} = DualArray{T, N}(A, zero.(A))
DualArray(A::AbstractArray{T, N}, B::AbstractArray{T, N}) where {T <: ReComp, N} = DualArray{T, N}(A, B)
DualArray(A::AbstractArray{T, N}) where {T <: ReComp, N} = DualArray{T, N}(A, zero.(A))

const DualVector{T} = DualArray{T, 1}
const DualMatrix{T} = DualArray{T, 2}

realpart(h::DualArray) = h.value
ɛpart(h::DualArray) = h.epsilon
realpart(x::AbstractArray) = x
ɛpart(x::AbstractArray) = zero.(typeof.(x))

isdualarray(::DualArray) = true
isdualarray(::AbstractArray) = false
isdualarray(::Dual) = true
isdualarray(::Number) = false

Base.size(A::DualArray) = size(realpart(A))

dualtype(::Type{Dual{T}}) where T = T

Base.similar(A::DualArray) = DualArray(similar(realpart(A)), similar(εpart(A)))
Base.similar(A::DualArray, dims::Vararg{Union{Integer, AbstractUnitRange}}) = DualArray(similar(realpart(A), dims), similar(εpart(A), dims...))
Base.similar(A::DualArray, dims::Tuple{Vararg{Union{Integer, AbstractUnitRange}}}) = DualArray(similar(realpart(A), dims), similar(εpart(A), dims))
Base.similar(A::AbstractArray, ::Type{T}, dims::Vararg{Union{Integer, AbstractUnitRange}}) where T <: Dual = DualArray(similar(realpart(A), dualtype(T), dims...), similar(εpart(A), dualtype(T), dims...))
# Base.similar(A::AbstractArray, ::Type{T}, dims::Tuple{Vararg{Union{Integer, AbstractUnitRange}}}) where T <: Dual = DualArray(similar(realpart(A), dualtype(T), dims...), similar(εpart(A), dualtype(T), dims))

Base.ones(::Type{T}, dims::Vararg{Union{Integer, AbstractUnitRange}}) where T <: Dual = DualArray(ones(dualtype(T), dims...))
Base.ones(::Type{T}, dims::Tuple{Vararg{Integer, N}}) where {T <: Dual, N} = DualArray(ones(dualtype(T), dims))
Base.ones(::Type{T}, dims::Tuple{Vararg{Base.OneTo, N}}) where {T <: Dual, N} = DualArray(ones(dualtype(T), dims))
Base.zeros(::Type{T}, dims::Vararg{Union{Integer, AbstractUnitRange}}) where T <: Dual = DualArray(zeros(dualtype(T), dims...))
Base.zeros(::Type{T}, dims::Tuple{Vararg{Integer, N}}) where {T <: Dual, N} = DualArray(zeros(dualtype(T), dims))
Base.zeros(::Type{T}, dims::Tuple{Vararg{Base.OneTo, N}}) where {T <: Dual, N} = DualArray(zeros(dualtype(T), dims))

Base.getindex(A::DualArray, inds...) = Dual(getindex(realpart(A), inds...), getindex(εpart(A), inds...))

function Base.setindex!(A::DualArray, X, inds...)
    setindex!(realpart(A), realpart(X), inds...)
    setindex!(εpart(A), εpart(X), inds...)
end

Base.:(==)(A::DualArray, B::DualArray) = realpart(A) == realpart(B)
Base.:(==)(A::DualArray, x::AbstractArray) = realpart(A) == x
Base.:(==)(x::AbstractArray, h::DualArray) = h == x

Base.isequal(A::DualArray, B::DualArray) = isequal(realpart(A), realpart(B)) && isequal(εpart(A), εpart(B))
Base.isequal(A::DualArray, x::AbstractArray) = isequal(realpart(A), x) && isequal(εpart(A), zero(x))
Base.isequal(x::AbstractArray, A::DualArray) = isequal(A, x)

Base.:-(A::DualArray) = DualArray(-realpart(A), -εpart(A))
Base.:-(A::DualArray, B::DualArray) = DualArray(realpart(A) - realpart(B), -εpart(B))
Base.:-(A::DualArray, B::DualArray) = DualArray(realpart(A) - realpart(B), εpart(A))
Base.:-(A::DualArray, B::DualArray) = DualArray(realpart(A) - realpart(B), εpart(A) - εpart(B))

Base.:+(A::DualArray, B::DualArray) = DualArray(realpart(A) + realpart(B), εpart(A) + εpart(B))

function apply_scalar(f, args::Vararg{Any, N}; kwargs...) where N
    if any(isdualarray, args)
        Dual(f(map(realpart, args)...; kwargs...), sum(map(i -> f(map(realpart, args[1 : i - 1])..., εpart(args[i]), map(realpart, args[i + 1 : end])...,; kwargs...), 1 : N)))
    else
        f(args...; kwargs...)
    end
end

function apply_linear(f, args::Vararg{Any, N}; kwargs...) where N
    if any(isdualarray, args)
        DualArray(f(map(realpart, args)...; kwargs...), sum(map(i -> f(map(realpart, args[1 : i - 1])..., εpart(args[i]), map(realpart, args[i + 1 : end])...,; kwargs...), 1 : N)))
    else
        f(args...; kwargs...)
    end
end

Base.:*(a::Dual, B::AbstractArray) = apply_linear(*, a, B)
Base.:*(A::AbstractArray, b::Dual) = apply_linear(*, A, b)
Base.:\(a::Dual, B::AbstractArray) = apply_linear(\, a, B)
Base.:/(A::AbstractArray, b::Dual) = apply_linear(/, A, b)

Base.:*(a::Number, B::DualArray) = apply_linear(*, a, B)
Base.:*(A::DualArray, b::Number) = apply_linear(*, A, b)
Base.:\(a::Number, B::DualArray) = apply_linear(\, a, B)
Base.:/(A::DualArray, b::Number) = apply_linear(/, A, b)

Base.:*(a::Dual, B::DualArray) = apply_linear(*, a, B)
Base.:*(A::DualArray, b::Dual) = apply_linear(*, A, b)
Base.:\(a::Dual, B::DualArray) = apply_linear(\, a, B)
Base.:/(A::DualArray, b::Dual) = apply_linear(/, A, b)

Base.:*(A::DualMatrix, B::AbstractVector) = apply_linear(*, A, B)
Base.:*(A::AbstractMatrix, B::DualVector) = apply_linear(*, A, B)
Base.:*(A::DualMatrix, B::DualVector) = apply_linear(*, A, B)

Base.:*(A::DualMatrix, B::AbstractMatrix) = apply_linear(*, A, B)
Base.:*(A::AbstractMatrix, B::DualMatrix) = apply_linear(*, A, B)
Base.:*(A::DualMatrix, B::DualMatrix) = apply_linear(*, A, B)

LinearAlgebra.dot(A::DualVector, B::AbstractVector) = apply_scalar(dot, A, B)
LinearAlgebra.dot(A::AbstractVector, B::DualVector) = apply_scalar(dot, A, B)
LinearAlgebra.dot(A::DualVector, B::DualVector) = apply_scalar(dot, A, B)

LinearAlgebra.dot(A::DualVector, B::AbstractMatrix, C::AbstractVector) = apply_scalar(dot, A, B, C)
LinearAlgebra.dot(A::AbstractVector, B::DualMatrix, C::AbstractVector) = apply_scalar(dot, A, B, C)
LinearAlgebra.dot(A::AbstractVector, B::AbstractMatrix, C::DualVector) = apply_scalar(dot, A, B, C)
LinearAlgebra.dot(A::DualVector, B::DualMatrix, C::AbstractVector) = apply_scalar(dot, A, B, C)
LinearAlgebra.dot(A::DualVector, B::AbstractMatrix, C::DualVector) = apply_scalar(dot, A, B, C)
LinearAlgebra.dot(A::AbstractVector, B::DualMatrix, C::DualVector) = apply_scalar(dot, A, B, C)
LinearAlgebra.dot(A::DualVector, B::DualMatrix, C::DualVector) = apply_scalar(dot, A, B, C)

Base.:*(A::Adjoint{T, <:AbstractVector} where T, B::DualMatrix) = dualein"i, ij -> j"(conj(parent(A)), B)
Base.:*(A::Adjoint{T, <:AbstractMatrix} where T, B::DualMatrix) = dualein"ij, ik -> jk"(conj(parent(A)), B)
Base.:*(A::DualMatrix, B::Adjoint{T, <:AbstractMatrix} where T) = dualein"ij, kj -> ik"(A, conj(parent(B)))

Base.:*(A::Adjoint{<: Number, <:AbstractVector}, B::DualVector) = dot(parent(A), B)
Base.:*(A::Adjoint{<: Number, <:DualVector}, B::AbstractVector) = dot(parent(A), B)
Base.:*(A::Adjoint{<: Number, <:DualVector}, B::DualVector) = dot(parent(A), B)

Base.:*(A::Adjoint{T, <:AbstractMatrix} where T, B::DualVector) = dualein"ij, i -> j"(conj(parent(A)), B)
Base.:*(A::Adjoint{T, <:DualMatrix} where T, B::AbstractVector) = dualein"ij, i -> j"(conj(parent(A)), B)
Base.:*(A::Adjoint{T, <:DualMatrix} where T, B::DualVector) = dualein"ij, i -> j"(conj(parent(A)), B)

LinearAlgebra.tr(A::DualMatrix) = Dual(tr(realpart(A)), tr(εpart(A)))

Base.conj(A::DualArray) = DualArray(conj(realpart(A)), conj(εpart(A)))
Base.reshape(A::DualArray, dims::Vararg{Int64, N}) where N = DualArray(reshape(realpart(A), dims), reshape(εpart(A), dims))
Base.reshape(A::DualArray, dims::Tuple{Vararg{Int64, N}}) where N = DualArray(reshape(realpart(A), dims), reshape(εpart(A), dims))