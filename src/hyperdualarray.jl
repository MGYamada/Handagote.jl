struct HyperDualArray{T <: ReComp, N} <: DenseArray{HyperDual{T}, N}
    value::Array{T, N}
    epsilon1::Array{T, N}
    epsilon2::Array{T, N}
    epsilon12::Array{T, N}
    function HyperDualArray{T, N}(A::AbstractArray{T, N}, B::AbstractArray{T, N}, C::AbstractArray{T, N}, D::AbstractArray{T, N}) where {T <: ReComp, N}
        @assert size(A) == size(B) == size(C) == size(D)
        new(convert(Array, A), convert(Array, B), convert(Array, C), convert(Array, D))
    end
end

HyperDualArray{T, N}(A::AbstractArray{Dual{T}, N}, B::AbstractArray{Dual{T}, N}) where {T <: ReComp, N} = HyperDualArray{T, N}(realpart(A), εpart(A), realpart(B), εpart(B))
HyperDualArray{T, N}(A::AbstractArray{Dual{T}, N}) where {T <: ReComp, N} = HyperDualArray{T, N}(A, DualArray(zero.(realpart(A))))
HyperDualArray(A::AbstractArray{Dual{T}, N}, B::AbstractArray{Dual{T}, N}) where {T <: ReComp, N} = HyperDualArray{T, N}(A, B)
HyperDualArray(A::AbstractArray{Dual{T}, N}) where {T <: ReComp, N} = HyperDualArray{T, N}(A, DualArray(zero.(realpart(A))))

HyperDualArray{T, N}(A::AbstractArray{T, N}) where {T <: ReComp, N} = HyperDualArray{T, N}(A, zero.(A), zero.(A), zero.(A))
HyperDualArray(A::AbstractArray{T, N}, B::AbstractArray{T, N}, C::AbstractArray{T, N}, D::AbstractArray{T, N}) where {T <: ReComp, N} = HyperDualArray{T, N}(A, B, C, D)
HyperDualArray(A::AbstractArray{T, N}) where {T <: ReComp, N} = HyperDualArray{T, N}(A, zero.(A), zero.(A), zero.(A))

const HyperDualVector{T} = HyperDualArray{T, 1}
const HyperDualMatrix{T} = HyperDualArray{T, 2}

realpart(h::HyperDualArray) = DualArray(h.value, h.epsilon1)
ɛpart(h::HyperDualArray) = DualArray(h.epsilon2, h.epsilon12)

hyperrealpart(h::HyperDualArray) = h.value
ɛ₁part(h::HyperDualArray) = h.epsilon1
ɛ₂part(h::HyperDualArray) = h.epsilon2
ɛ₁ε₂part(h::HyperDualArray) = h.epsilon12

hyperrealpart(d::DualArray) = realpart(d)
ɛ₁part(d::DualArray) = εpart(d)
ɛ₂part(d::DualArray) = zero.(typeof.(realpart(d)))
ɛ₁ε₂part(d::DualArray) = zero.(typeof.(realpart(d)))

hyperrealpart(x::AbstractArray) = x
ɛ₁part(x::AbstractArray) = zero.(typeof.(x))
ɛ₂part(x::AbstractArray) = zero.(typeof.(x))
ɛ₁ε₂part(x::AbstractArray) = zero.(typeof.(x))

ishyperdualarray(::HyperDualArray) = true
ishyperdualarray(::AbstractArray) = false
isdualarray(::HyperDual) = true
isdualarray(::Number) = false

Base.size(A::HyperDualArray) = size(hyperrealpart(A))

hyperdualtype(::Type{HyperDual{T}}) where T = T

Base.similar(A::HyperDualArray) = HyperDualArray(similar(realpart(A)), similar(εpart(A)))
Base.similar(A::HyperDualArray, dims::Vararg{Union{Integer, AbstractUnitRange}}) = HyperDualArray(similar(hyperrealpart(A), dims))
Base.similar(A::HyperDualArray, dims::Tuple{Vararg{Union{Integer, AbstractUnitRange}}}) = HyperDualArray(similar(hyperrealpart(A), dims))
Base.similar(A::AbstractArray, ::Type{HyperDual{T}}, dims::Vararg{Union{Integer, AbstractUnitRange}}) where T = HyperDualArray(similar(hyperrealpart(A), T, dims...))
# Base.similar(A::AbstractArray, ::Type{HyperDual{T}}, dims::Tuple{Vararg{Union{Integer, AbstractUnitRange}}}) where T = HyperDualArray(similar(hyperrealpart(A), T, dims...))

Base.ones(::Type{T}, dims::Vararg{Union{Integer, AbstractUnitRange}}) where T <: HyperDual = DualArray(ones(hyperdualtype(T), dims...))
Base.ones(::Type{T}, dims::Tuple{Vararg{Integer, N}}) where {T <: HyperDual, N} = DualArray(ones(hyperdualtype(T), dims))
Base.ones(::Type{T}, dims::Tuple{Vararg{Base.OneTo, N}}) where {T <: HyperDual, N} = DualArray(ones(hyperdualtype(T), dims))
Base.zeros(::Type{T}, dims::Vararg{Union{Integer, AbstractUnitRange}}) where T <: HyperDual = DualArray(zeros(hyperdualtype(T), dims...))
Base.zeros(::Type{T}, dims::Tuple{Vararg{Integer, N}}) where {T <: HyperDual, N} = DualArray(zeros(hyperdualtype(T), dims))
Base.zeros(::Type{T}, dims::Tuple{Vararg{Base.OneTo, N}}) where {T <: HyperDual, N} = DualArray(zeros(hyperdualtype(T), dims))

Base.getindex(A::HyperDualArray, inds...) = HyperDual(getindex(hyperrealpart(A), inds...), getindex(ɛ₁part(A), inds...), getindex(ɛ₂part(A), inds...), getindex(ɛ₁ε₂part(A), inds...))

function Base.setindex!(A::HyperDualArray, X, inds...)
    setindex!(hyperrealpart(A), hyperrealpart(X), inds...)
    setindex!(ɛ₁part(A), ɛ₁part(X), inds...)
    setindex!(ɛ₂part(A), ɛ₂part(X), inds...)
    setindex!(ɛ₁ε₂part(A), ɛ₁ε₂part(X), inds...)
end

Base.:(==)(A::HyperDualArray, B::HyperDualArray) = realpart(A) == realpart(B)
Base.:(==)(A::HyperDualArray, x::AbstractArray) = realpart(A) == x
Base.:(==)(x::AbstractArray, h::HyperDualArray) = h == x

Base.isequal(A::HyperDualArray, B::HyperDualArray) = isequal(realpart(A), realpart(B)) && isequal(εpart(A), εpart(B))
Base.isequal(A::HyperDualArray, x::AbstractArray) = isequal(realpart(A), x) && isequal(εpart(A), zero(x))
Base.isequal(x::AbstractArray, A::HyperDualArray) = isequal(A, x)

Base.:-(A::HyperDualArray) = HyperDualArray(-realpart(A), -εpart(A))
Base.:-(A::HyperDualArray, B::HyperDualArray) = HyperDualArray(realpart(A) - realpart(B), εpart(A) - εpart(B))

Base.:+(A::HyperDualArray, B::HyperDualArray) = HyperDualArray(realpart(A) + realpart(B), εpart(A) + εpart(B))

function apply_scalar_hyper(f, args::Vararg{Any, N}; kwargs...) where N
    if any(ishyperdualarray, args)
        HyperDual(f(map(hyperrealpart, args)...; kwargs...), sum(i -> f(map(hyperrealpart, args[1 : i - 1])..., ɛ₁part(args[i]), map(hyperrealpart, args[i + 1 : end])...,; kwargs...), 1 : N),
        sum(i -> f(map(hyperrealpart, args[1 : i - 1])..., ɛ₂part(args[i]), map(hyperrealpart, args[i + 1 : end])...,; kwargs...), 1 : N),
        sum(i -> sum(j -> i < j ? f(map(hyperrealpart, args[1 : i - 1])..., ɛ₁part(args[i]), map(hyperrealpart, args[i + 1 : j - 1])..., ɛ₂part(args[j]), map(hyperrealpart, args[j + 1 : end])...,; kwargs...) :
        (i > j ? f(map(hyperrealpart, args[1 : j - 1])..., ɛ₂part(args[j]), map(hyperrealpart, args[j + 1 : i - 1])..., ɛ₁part(args[i]), map(hyperrealpart, args[i + 1 : end])...,; kwargs...) :
        f(map(hyperrealpart, args[1 : i - 1])..., ɛ₁ɛ₂part(args[i]), map(hyperrealpart, args[i + 1 : end])...,; kwargs...)), 1 : N), 1 : N))
    else
        apply_scalar(f, args...; kwargs...) 
    end
end

function apply_linear_hyper(f, args::Vararg{Any, N}; kwargs...) where N
    if any(ishyperdualarray, args)
        HyperDualArray(f(map(hyperrealpart, args)...; kwargs...), sum(i -> f(map(hyperrealpart, args[1 : i - 1])..., ɛ₁part(args[i]), map(hyperrealpart, args[i + 1 : end])...,; kwargs...), 1 : N),
        sum(i -> f(map(hyperrealpart, args[1 : i - 1])..., ɛ₂part(args[i]), map(hyperrealpart, args[i + 1 : end])...,; kwargs...), 1 : N),
        sum(i -> sum(j -> i < j ? f(map(hyperrealpart, args[1 : i - 1])..., ɛ₁part(args[i]), map(hyperrealpart, args[i + 1 : j - 1])..., ɛ₂part(args[j]), map(hyperrealpart, args[j + 1 : end])...,; kwargs...) :
        (i > j ? f(map(hyperrealpart, args[1 : j - 1])..., ɛ₂part(args[j]), map(hyperrealpart, args[j + 1 : i - 1])..., ɛ₁part(args[i]), map(hyperrealpart, args[i + 1 : end])...,; kwargs...) :
        f(map(hyperrealpart, args[1 : i - 1])..., ɛ₁ɛ₂part(args[i]), map(hyperrealpart, args[i + 1 : end])...,; kwargs...)), 1 : N), 1 : N))
    else
        apply_linear(f, args...; kwargs...) 
    end
end

Base.:*(a::HyperDual, B::AbstractArray) = apply_linear_hyper(*, a, B)
Base.:*(A::AbstractArray, b::HyperDual) = apply_linear_hyper(*, A, b)
Base.:\(a::HyperDual, B::AbstractArray) = apply_linear_hyper(\, a, B)
Base.:/(A::AbstractArray, b::HyperDual) = apply_linear_hyper(/, A, b)

Base.:*(a::Number, B::HyperDualArray) = apply_linear_hyper(*, a, B)
Base.:*(A::HyperDualArray, b::Number) = apply_linear_hyper(*, A, b)
Base.:\(a::Number, B::HyperDualArray) = apply_linear_hyper(\, a, B)
Base.:/(A::HyperDualArray, b::Number) = apply_linear_hyper(/, A, b)

Base.:*(a::HyperDual, B::HyperDualArray) = apply_linear_hyper(*, a, B)
Base.:*(A::HyperDualArray, b::HyperDual) = apply_linear_hyper(*, A, b)
Base.:\(a::HyperDual, B::HyperDualArray) = apply_linear_hyper(\, a, B)
Base.:/(A::HyperDualArray, b::HyperDual) = apply_linear_hyper(/, A, b)

Base.:*(A::HyperDualMatrix, B::AbstractVector) = apply_linear_hyper(*, A, B)
Base.:*(A::AbstractMatrix, B::HyperDualVector) = apply_linear_hyper(*, A, B)
Base.:*(A::HyperDualMatrix, B::HyperDualVector) = apply_linear_hyper(*, A, B)

Base.:*(A::HyperDualMatrix, B::AbstractMatrix) = apply_linear_hyper(*, A, B)
Base.:*(A::AbstractMatrix, B::HyperDualMatrix) = apply_linear_hyper(*, A, B)
Base.:*(A::HyperDualMatrix, B::HyperDualMatrix) = apply_linear_hyper(*, A, B)

Base.:*(A::Adjoint{T, <:AbstractVector} where T, B::HyperDualMatrix) = hyperdualein"i, ij -> j"(conj(parent(A)), B)
Base.:*(A::Adjoint{T, <:AbstractMatrix} where T, B::HyperDualMatrix) = hyperdualein"ij, ik -> jk"(conj(parent(A)), B)
Base.:*(A::HyperDualMatrix, B::Adjoint{T, <:AbstractMatrix} where T) = hyperdualein"ij, kj -> ik"(A, conj(parent(B)))

Base.:*(A::Adjoint{T, <:AbstractMatrix} where T, B::HyperDualVector) = hyperdualein"ij, i -> j"(conj(parent(A)), B)
Base.:*(A::Adjoint{T, <:HyperDualMatrix} where T, B::AbstractVector) = hyperdualein"ij, i -> j"(conj(parent(A)), B)
Base.:*(A::Adjoint{T, <:HyperDualMatrix} where T, B::HyperDualVector) = hyperdualein"ij, i -> j"(conj(parent(A)), B)

LinearAlgebra.tr(A::HyperDualMatrix) = HyperDual(tr(hyperrealpart(A)), tr(ɛ₁part(A)), tr(ɛ₂part(A)), tr(ɛ₁ε₂part(A)))

Base.conj(A::HyperDualArray) = HyperDualArray(conj(hyperrealpart(A)), conj(ɛ₁part(A)), conj(ɛ₂part(A)), conj(ɛ₁ε₂part(A)))
Base.reshape(A::HyperDualArray, dims::Vararg{Int64, N}) where N = HyperDualArray(reshape(realpart(A), dims), reshape(εpart(A), dims))
Base.reshape(A::HyperDualArray, dims::Tuple{Vararg{Int64, N}}) where N = HyperDualArray(reshape(realpart(A), dims), reshape(εpart(A), dims))