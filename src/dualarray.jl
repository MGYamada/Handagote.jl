abstract type AbstractDualArray{T, N} end

const AbstractDualVector{T} = AbstractDualArray{T, 1}
const AbstractDualMatrix{T} = AbstractDualArray{T, 2}

mutable struct DualArray{T <: NumberOrDual, N} <: AbstractDualArray{Dual{T}, N}
    value::Array{T, N}
    epsilon::Array{T, N}
    function DualArray{T, N}(A::Array{T, N}, B::Array{T, N}) where {T <: NumberOrDual, N}
        @assert size(A) == size(B)
        new(A, B)
    end
end

const AbstractArrayOrDualArray = Union{AbstractArray, AbstractDualArray}

DualArray{T, N}(A::Array{T, N}) where {T <: NumberOrDual, N} = DualArray{T, N}(A, zero.(A))
DualArray(A::Array{T, N}, B::Array{T, N}) where {T <: NumberOrDual, N} = DualArray{T, N}(A, B)
DualArray(A::Array{T, N}) where {T <: NumberOrDual, N} = DualArray{T, N}(A, zero.(A))

const DualVector{T} = DualArray{T, 1}
const DualMatrix{T} = DualArray{T, 2}

realpart(h::AbstractDualArray) = h.value
ɛpart(h::AbstractDualArray) = h.epsilon
realpart(x::AbstractArrayOrDualArray) = x
ɛpart(x::AbstractArrayOrDualArray) = zero.(typeof.(x))

isdualarray(::AbstractDualArray) = true
isdualarray(::AbstractArrayOrDualArray) = false

Base.size(A::AbstractDualArray) = size(realpart(A))
Base.size(A::AbstractDualArray, i::Int) = size(realpart(A), i)

dualtype(::Type{Dual{T}}) where T = T

Base.copy(A::AbstractDualArray) = DualArray(copy(realpart(A)), copy(εpart(A)))

Base.similar(A::AbstractDualArray) = DualArray(similar(realpart(A)), similar(εpart(A)))
Base.similar(A::AbstractDualArray, dims::Vararg{Union{Integer, AbstractUnitRange}}) = DualArray(similar(realpart(A), dims), similar(εpart(A), dims...))
Base.similar(A::AbstractDualArray, dims::Tuple{Vararg{Union{Integer, AbstractUnitRange}}}) = DualArray(similar(realpart(A), dims), similar(εpart(A), dims))
Base.similar(A::AbstractArrayOrDualArray, ::Type{T}, dims::Vararg{Union{Integer, AbstractUnitRange}}) where T <: AbstractDual = DualArray(similar(realpart(A), dualtype(T), dims...), similar(εpart(A), dualtype(T), dims...))
Base.similar(A::AbstractArrayOrDualArray, ::Type{T}, dims::Tuple{Vararg{Union{Integer, AbstractUnitRange}}}) where T <: AbstractDual = DualArray(similar(realpart(A), dualtype(T), dims...), similar(εpart(A), dualtype(T), dims))

Base.ones(::Type{T}, dims::Vararg{Union{Integer, AbstractUnitRange}}) where T <: AbstractDual = DualArray(ones(dualtype(T), dims...))
Base.ones(::Type{T}, dims::Tuple{Vararg{Integer, N}}) where {T <: AbstractDual, N} = DualArray(ones(dualtype(T), dims))
Base.ones(::Type{T}, dims::Tuple{Vararg{Base.OneTo, N}}) where {T <: AbstractDual, N} = DualArray(ones(dualtype(T), dims))
Base.zeros(::Type{T}, dims::Vararg{Union{Integer, AbstractUnitRange}}) where T <: AbstractDual = DualArray(zeros(dualtype(T), dims...))
Base.zeros(::Type{T}, dims::Tuple{Vararg{Integer, N}}) where {T <: AbstractDual, N} = DualArray(zeros(dualtype(T), dims))
Base.zeros(::Type{T}, dims::Tuple{Vararg{Base.OneTo, N}}) where {T <: AbstractDual, N} = DualArray(zeros(dualtype(T), dims))

Base.getindex(A::AbstractDualArray, inds...) = Dual(getindex(realpart(A), inds...), getindex(εpart(A), inds...))

function Base.setindex!(A::AbstractDualArray, X, inds...)
    setindex!(realpart(A), realpart(X), inds...)
    setindex!(εpart(A), εpart(X), inds...)
end

Base.:(==)(A::AbstractDualArray, B::AbstractDualArray) = realpart(A) == realpart(B)
Base.:(==)(A::AbstractDualArray, x::AbstractArrayOrDualArray) = realpart(A) == x
Base.:(==)(x::AbstractArrayOrDualArray, h::AbstractDualArray) = h == x

Base.isequal(A::AbstractDualArray, B::AbstractDualArray) = isequal(realpart(A), realpart(B)) && isequal(εpart(A), εpart(B))
Base.isequal(A::AbstractDualArray, x::AbstractArrayOrDualArray) = isequal(realpart(A), x) && isequal(εpart(A), zero(x))
Base.isequal(x::AbstractArrayOrDualArray, A::AbstractDualArray) = isequal(A, x)

Base.:-(A::AbstractDualArray) = DualArray(-realpart(A), -εpart(A))
Base.:-(A::AbstractDualArray, B::AbstractDualArray) = DualArray(realpart(A) - realpart(B), εpart(A) - εpart(B))

Base.:+(A::AbstractDualArray, B::AbstractDualArray) = DualArray(realpart(A) + realpart(B), εpart(A) + εpart(B))

function apply_linear(f, args::Vararg{Any, N}; kwargs...) where N
    DualArray(f(map(realpart, args)...; kwargs...), sum(map(i -> f(map(realpart, args[1 : i - 1])..., εpart(args[i]), map(realpart, args[i + 1 : end])...,; kwargs...), 1 : N)))
end

Base.:*(a::Number, B::AbstractDualArray) = apply_linear(*, a, B)
Base.:*(A::AbstractDualArray, b::Number) = apply_linear(*, A, b)
Base.:\(a::Number, B::AbstractDualArray) = apply_linear(\, a, B)
Base.:/(A::AbstractDualArray, b::Number) = apply_linear(/, A, b)

Base.:*(A::AbstractDualMatrix, B::AbstractMatrix) = apply_linear(*, A, B)
Base.:*(A::AbstractMatrix, B::AbstractDualMatrix) = apply_linear(*, A, B)
Base.:*(A::AbstractDualMatrix, B::AbstractDualMatrix) = apply_linear(*, A, B)

LinearAlgebra.tr(A::AbstractDualMatrix) = Dual(tr(realpart(A)), tr(εpart(A)))

Base.conj(A::AbstractDualArray) = DualArray(conj(realpart(A)), conj(εpart(A)))
Base.reshape(A::AbstractDualArray, dims::Vararg{Int64, N}) where N = DualArray(reshape(realpart(A), dims), reshape(εpart(A), dims))
Base.reshape(A::AbstractDualArray, dims::Tuple{Vararg{Int64, N}}) where N = DualArray(reshape(realpart(A), dims), reshape(εpart(A), dims))