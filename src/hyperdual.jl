mutable struct HyperDual{T <: ReComp} <: Number
    value::T
    epsilon1::T
    epsilon2::T
    epsilon12::T
end

HyperDual(x::Dual, y::Dual) = HyperDual(realpart(x), εpart(x), realpart(y), εpart(y))
HyperDual{T}(x::Dual, y::Dual) where T <: ReComp = HyperDual{T}(T(realpart(x)), T(εpart(x)), T(realpart(y)), T(εpart(y)))
HyperDual(x::Dual) = HyperDual(realpart(x), εpart(x), zero(x), zero(x))
HyperDual{T}(x::Dual) where T <: ReComp = HyperDual{T}(T(realpart(x)), T(εpart(x)), zero(T), zero(T))

HyperDual(x::S, y::T, z::U, w::V) where {S <: ReComp, T <: ReComp, U <: ReComp, V <: ReComp} = HyperDual(promote(x, y, z, w)...)
HyperDual(x::ReComp) = HyperDual(x, zero(x), zero(x), zero(x))
HyperDual{T}(x::ReComp) where T <: ReComp = HyperDual{T}(T(x), zero(T), zero(T), zero(T))

HyperDual() = HyperDual(false, false, false, false)
const ɛ₁ = HyperDual(false, true, false, false)
const ɛ₂ = HyperDual(false, false, true, false)
const ε₁ɛ₂ = HyperDual(false, false, false, true)

Base.convert(::Type{HyperDual{T}}, h::HyperDual{T}) where {T<:ReComp} = h
Base.convert(::Type{HyperDual{T}}, h::HyperDual) where {T<:ReComp} = HyperDual{T}(convert(T, hyperrealpart(h)), convert(T, ɛ₁part(h)), convert(T, ɛ₂part(h)), convert(T, ɛ₁ε₂part(h)))
Base.convert(::Type{HyperDual{T}}, x::Number) where {T<:ReComp} = HyperDual{T}(convert(T, x), convert(T, 0), convert(T, 0), convert(T, 0))

Base.convert(::Type{HyperDual{T}}, h::Dual{T}) where {T<:ReComp} = HyperDual(h)
Base.convert(::Type{HyperDual{T}}, h::Dual) where {T<:ReComp} = HyperDual{T}(convert(T, hyperrealpart(h)), convert(T, ɛ₁part(h)))

Base.convert(::Type{T}, h::HyperDual) where {T <: Dual} = (ɛ₂part(h) == 0 && ɛ₁ε₂part(h) == 0 ? convert(T, realpart(value(h))) : throw(InexactError()))
Base.convert(::Type{T}, h::HyperDual) where {T <: ReComp} = (ɛ₁part(h) == 0 && ɛ₂part(h) == 0 && ɛ₁ε₂part(h) == 0 ? convert(T, hyperrealpart(h)) : throw(InexactError()))

Base.promote_rule(::Type{HyperDual{T}}, ::Type{HyperDual{S}}) where {T<:ReComp, S<:ReComp} = HyperDual{promote_type(T, S)}
Base.promote_rule(::Type{HyperDual{T}}, ::Type{S}) where {T<:ReComp, S<:ReComp} = HyperDual{promote_type(T, S)}
Base.promote_rule(::Type{HyperDual{T}}, ::Type{T}) where {T<:ReComp} = HyperDual{T}

Base.widen(::Type{HyperDual{T}}) where T = HyperDual{widen(T)}

realpart(h::HyperDual) = Dual(h.value, h.epsilon1)
ɛpart(h::HyperDual) = Dual(h.epsilon2, h.epsilon12)

hyperrealpart(h::HyperDual) = h.value
ɛ₁part(h::HyperDual) = h.epsilon1
ɛ₂part(h::HyperDual) = h.epsilon2
ɛ₁ε₂part(h::HyperDual) = h.epsilon12

hyperrealpart(d::Dual) = realpart(d)
ɛ₁part(d::Dual) = ɛpart(d)
ɛ₂part(d::Dual) = zero(typeof(realpart(d)))
ɛ₁ε₂part(x::Dual) = zero(typeof(realpart(d)))

hyperrealpart(x::Number) = x
ɛ₁part(x::Number) = zero(typeof(x))
ɛ₂part(x::Number) = zero(typeof(x))
ɛ₁ε₂part(x::Number) = zero(typeof(x))

Base.isnan(h::HyperDual) = isnan(hyperrealpart(h))
Base.isinf(h::HyperDual) = isinf(hyperrealpart(h))
Base.isfinite(h::HyperDual) = isfinite(hyperrealpart(h))
ishyperdual(x::HyperDual) = true
ishyperdual(x::Number) = false
Base.eps(h::HyperDual) = eps(hyperrealpart(h))
Base.eps(::Type{HyperDual{T}}) where T = eps(T)