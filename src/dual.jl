abstract type AbstractDual{T} end

const NumberOrDual = Union{Number, AbstractDual}

mutable struct Dual{T <: NumberOrDual} <: AbstractDual{T}
    value::T
    epsilon::T
end

Dual(x::S, y::T) where {S <: NumberOrDual, T <: NumberOrDual} = Dual(promote(x, y)...)
Dual(x::NumberOrDual) = Dual(x, zero(x))
Dual{T}(x::NumberOrDual) where T <: NumberOrDual = Dual{T}(T(x), zero(T))
Dual() = Dual(false, false)

const ɛ = Dual(false, true)

Base.convert(::Type{AbstractDual{T}}, h::AbstractDual{T}) where {T <: NumberOrDual} = h
Base.convert(::Type{AbstractDual{T}}, h::AbstractDual) where {T <: NumberOrDual} = Dual{T}(convert(T, realpart(h)), convert(T, ɛpart(h)))
Base.convert(::Type{AbstractDual{T}}, x::NumberOrDual) where {T <: NumberOrDual} = Dual{T}(convert(T, x), convert(T, 0))
Base.convert(::Type{T}, h::AbstractDual) where {T <: NumberOrDual} = (ɛpart(h) == 0 ? convert(T, realpart(h)) : throw(InexactError()))

Base.promote_rule(::Type{AbstractDual{T}}, ::Type{Dual{S}}) where {T <: NumberOrDual, S <: NumberOrDual} = Dual{promote_type(T, S)}
Base.promote_rule(::Type{AbstractDual{T}}, ::Type{S}) where {T <: NumberOrDual, S <: NumberOrDual} = Dual{promote_type(T, S)}
Base.promote_rule(::Type{AbstractDual{T}}, ::Type{T}) where {T <: NumberOrDual} = Dual{T}

Base.widen(::Type{AbstractDual{T}}) where T = Dual{widen(T)}

realpart(h::AbstractDual) = h.value
ɛpart(h::AbstractDual) = h.epsilon
realpart(x::NumberOrDual) = x
ɛpart(x::NumberOrDual) = zero(typeof(x))

Base.isnan(h::AbstractDual) = isnan(realpart(h))
Base.isinf(h::AbstractDual) = isinf(realpart(h))
Base.isfinite(h::AbstractDual) = isfinite(realpart(h))
isdual(x::AbstractDual) = true
isdual(x::NumberOrDual) = false
Base.eps(h::AbstractDual) = eps(realpart(h))
Base.eps(::Type{AbstractDual{T}}) where T = eps(T)

Base.convert(::Type{Dual}, h::AbstractDual) = h
Base.convert(::Type{Dual}, x::NumberOrDual) = Dual(x)

Base.:(==)(h₁::AbstractDual, h₂::AbstractDual) = realpart(h₁) == realpart(h₂)
Base.:(==)(h::AbstractDual, x::NumberOrDual) = realpart(h) == x
Base.:(==)(x::NumberOrDual, h::AbstractDual) = h == x

Base.isequal(h₁::AbstractDual, h₂::AbstractDual) = isequal(realpart(h₁), realpart(h₂)) && isequal(εpart(h₁), εpart(h₂))
Base.isequal(h::AbstractDual, x::NumberOrDual) = isequal(realpart(h), x) && isequal(εpart(h), zero(x))
Base.isequal(x::NumberOrDual, h::AbstractDual) = isequal(h, x)

Base.isless(h₁::AbstractDual{T}, h₂::AbstractDual{T}) where {T <: Real} = realpart(h₁) < realpart(h₂)
Base.isless(h₁::Real, h₂::AbstractDual{<:Real}) = h₁ < realpart(h₂)
Base.isless(h₁::AbstractDual{<:Real}, h₂::Real) = realpart(h₁) < h₂

function Base.hash(h::AbstractDual) # Maybe this works
    x = hash(realpart(h))
    if isequal(h, realpart(h))
        return x
    else
        y = hash(εpart(h))
        return hash(x, hash(y))
    end
end

Base.float(h::Union{AbstractDual{T}, AbstractDual{Complex{T}}}) where {T<:AbstractFloat} = h
Base.complex(h::AbstractDual{<:Complex}) = h
Base.complex(::Type{AbstractDual{T}}) where {T} = Dual{complex(T)}

Base.floor(h::AbstractDual) = floor(realpart(h))
Base.ceil(h::AbstractDual)  = ceil(realpart(h))
Base.trunc(h::AbstractDual) = trunc(realpart(h))
Base.round(h::AbstractDual) = round(realpart(h))
Base.floor(::Type{T}, h::AbstractDual) where {T<:Real} = floor(T, realpart(h))
Base.ceil( ::Type{T}, h::AbstractDual) where {T<:Real} =  ceil(T, realpart(h))
Base.trunc(::Type{T}, h::AbstractDual) where {T<:Real} = trunc(T, realpart(h))
Base.round(::Type{T}, h::AbstractDual) where {T<:Real} = round(T, realpart(h))

for op in (:real, :imag, :conj, :float, :complex)
    @eval Base.$op(z::AbstractDual) = Dual($op(realpart(z)), $op(εpart(z)))
end

Base.abs(z::AbstractDual) = sqrt(abs2(z))
Base.abs2(z::AbstractDual) = real(conj(z) * z)

Base.real(z::AbstractDual{<:Real}) = z
Base.abs(z::AbstractDual{<:Real}) = z ≥ 0 ? z : -z

Base.angle(z::AbstractDual{<:Real}) = z ≥ 0 ? zero(z) : one(z) * π
function Base.angle(z::AbstractDual{Complex{T}}) where T <: Real
    if z == 0
        if imag(εpart(z)) == 0
            Dual(zero(T), zero(T))
        else
            Dual(zero(T), convert(T, Inf))
        end
    else
        real(log(sign(z)) / im)
    end
end

Base.flipsign(x::AbstractDual, y::AbstractDual) = y == 0 ? flipsign(x, εpart(y)) : flipsign(x, realpart(y))
Base.flipsign(x, y::AbstractDual) = y == 0 ? flipsign(x, εpart(y)) : flipsign(x, realpart(y))
Base.flipsign(x::AbstractDual, y) = dual(flipsign(realpart(x), y), flipsign(εpart(x), y))

Base.:+(h₁::AbstractDual, h₂::AbstractDual) = Dual(realpart(h₁) + realpart(h₂), εpart(h₁) + εpart(h₂))
Base.:+(n::NumberOrDual, h::AbstractDual) = Dual(n + realpart(h), εpart(h))
Base.:+(h::AbstractDual, n::NumberOrDual) = n + h

Base.:-(h::AbstractDual) = Dual(-realpart(h), -εpart(h))
Base.:-(h₁::AbstractDual, h₂::AbstractDual) = Dual(realpart(h₁) - realpart(h₂), εpart(h₁) - εpart(h₂))
Base.:-(n::NumberOrDual, h::AbstractDual) = Dual(n - realpart(h), -εpart(h))
Base.:-(h::AbstractDual, n::NumberOrDual) = Dual(realpart(h) - n, εpart(h))

# avoid ambiguous definition
Base.:*(x::Bool, h::AbstractDual) = ifelse(x, h, ifelse(signbit(real(realpart(h))) == 0, zero(h), -zero(h)))
Base.:*(h::AbstractDual, x::Bool) = x * h

function Base.:*(h₁::AbstractDual, h₂::AbstractDual)
    x, y = realpart(h₁), εpart(h₁)
    a, b = realpart(h₂), εpart(h₂)
    return Dual(a*x, muladd(a, y, b*x))
end
Base.:*(n::NumberOrDual, h::AbstractDual) = Dual(n*realpart(h), n*εpart(h))
Base.:*(h::AbstractDual, n::NumberOrDual) = n * h

Base.one(h::AbstractDual) = Dual(one(realpart(h)))

@inline Base.literal_pow(::typeof(^), x::AbstractDual, ::Val{0}) = one(typeof(x))
@inline Base.literal_pow(::typeof(^), x::AbstractDual, ::Val{1}) = x
@inline Base.literal_pow(::typeof(^), x::AbstractDual, ::Val{2}) = x * x
@inline Base.literal_pow(::typeof(^), x::AbstractDual, ::Val{3}) = x * x * x

function Base.:/(h₁::AbstractDual, h₂::AbstractDual)
    x, y = realpart(h₁), εpart(h₁)
    a, b = realpart(h₂), εpart(h₂)
    return Dual(x / a, y / a - b * x / a ^ 2)
end
function Base.:/(n::NumberOrDual, h::AbstractDual)
    x, y = realpart(h), ε₁part(h)
    return Dual(n / x, -n * y / x ^ 2)
end
Base.:/(h::AbstractDual, n::NumberOrDual) = Dual(realpart(h) / n, εpart(h) / n)

Base.mod(h::AbstractDual, n::NumberOrDual) = Dual(mod(realpart(h), n), εpart(h))

# power

function to_nanmath(x::Expr)
    if x.head == :call
        funsym = Expr(:.,:NaNMath,Base.Meta.quot(x.args[1]))
        return Expr(:call,funsym,[to_nanmath(z) for z in x.args[2:end]]...)
    else
        return Expr(:call,[to_nanmath(z) for z in x.args]...)
    end
end
to_nanmath(x) = x

Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, i::AbstractDual) = checkindex(Bool, inds, realpart(h))