mutable struct HyperDual{T <: NumberOrDual} <: AbstractDual{Dual{T}}
    value::T
    epsilon1::T
    epsilon2::T
    epsilon12::T
end

HyperDual(x::S, y::T, z::U, w::V) where {S <: NumberOrDual, T <: NumberOrDual, U <: NumberOrDual, V <: NumberOrDual} = HyperDual(promote(x, y, z, w)...)
HyperDual(x::Dual) = HyperDual(realpart(x), εpart(x), zero(x), zero(x))
HyperDual{T}(x::Dual) where T <: NumberOrDual = HyperDual{T}(T(realpart(x)), T(εpart(x)), zero(T), zero(T))
HyperDual(x::NumberOrDual) = HyperDual(x, zero(x), zero(x), zero(x))
HyperDual{T}(x::NumberOrDual) where T <: NumberOrDual = HyperDual{T}(T(x), zero(T), zero(T), zero(T))

HyperDual() = HyperDual(false, false, false, false)
const ɛ₁ = HyperDual(false, true, false, false)
const ɛ₂ = HyperDual(false, false, true, false)
const ε₁ɛ₂ = HyperDual(false, false, false, true)

realpart(h::HyperDual) = Dual(h.value, h.epsilon1)
ɛpart(h::HyperDual) = Dual(h.epsilon2, h.epsilon12)

hyperrealpart(h::HyperDual) = h.value
ɛ₁part(h::HyperDual) = h.epsilon1
ɛ₂part(h::HyperDual) = h.epsilon2
ɛ₁ε₂part(h::HyperDual) = h.epsilon12

hyperrealpart(d::Dual) = realpart(d)
ɛ₁part(d::Dual) = εpart(d)
ɛ₂part(d::Dual) = zero(typeof(realpart(d)))
ɛ₁ε₂part(x::Dual) = zero(typeof(realpart(d)))

hyperrealpart(x::NumberOrDual) = x
ɛ₁part(x::NumberOrDual) = zero(typeof(x))
ɛ₂part(x::NumberOrDual) = zero(typeof(x))
ɛ₁ε₂part(x::NumberOrDual) = zero(typeof(x))

Base.isnan(h::HyperDual) = isnan(hyperrealpart(h))
Base.isinf(h::HyperDual) = isinf(hyperrealpart(h))
Base.isfinite(h::HyperDual) = isfinite(hyperrealpart(h))
ishyperdual(x::HyperDual) = true
ishyperdual(x::NumberOrDual) = false
Base.eps(h::HyperDual) = eps(hyperrealpart(h))
Base.eps(::Type{HyperDual{T}}) where T = eps(T)