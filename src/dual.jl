const ReComp = Union{Real, Complex}

struct Dual{T <: ReComp} <: Number
    value::T
    epsilon::T
end

Dual(x::S, y::T) where {S <: ReComp, T <: ReComp} = Dual(promote(x, y)...)
Dual(x::ReComp) = Dual(x, zero(x))
Dual{T}(x::ReComp) where T <: ReComp = Dual{T}(T(x), zero(T))
Dual() = Dual(false, false)

const ɛ = Dual(false, true)

Base.convert(::Type{Dual{T}}, h::Dual{T}) where {T <: ReComp} = h
Base.convert(::Type{Dual{T}}, h::Dual) where {T <: ReComp} = Dual{T}(convert(T, realpart(h)), convert(T, ɛpart(h)))
Base.convert(::Type{Dual{T}}, x::ReComp) where {T <: ReComp} = Dual{T}(convert(T, x), convert(T, 0))
Base.convert(::Type{T}, h::Dual) where {T <: ReComp} = (ɛpart(h) == 0 ? convert(T, realpart(h)) : throw(InexactError()))

Base.promote_rule(::Type{Dual{T}}, ::Type{Dual{S}}) where {T <: ReComp, S <: ReComp} = Dual{promote_type(T, S)}
Base.promote_rule(::Type{Dual{T}}, ::Type{S}) where {T <: ReComp, S <: ReComp} = Dual{promote_type(T, S)}
Base.promote_rule(::Type{Dual{T}}, ::Type{T}) where {T <: ReComp} = Dual{T}

Base.widen(::Type{Dual{T}}) where T = Dual{widen(T)}

realpart(h::Dual) = h.value
ɛpart(h::Dual) = h.epsilon
realpart(x::Number) = x
ɛpart(x::Number) = zero(typeof(x))

Base.isnan(h::Dual) = isnan(realpart(h))
Base.isinf(h::Dual) = isinf(realpart(h))
Base.isfinite(h::Dual) = isfinite(realpart(h))
isdual(x::Dual) = true
isdual(x::Number) = false
Base.eps(h::Dual) = eps(realpart(h))
Base.eps(::Type{Dual{T}}) where T = eps(T)

Base.convert(::Type{Dual}, h::Dual) = h
Base.convert(::Type{Dual}, x::Number) = Dual(x)

Base.:(==)(h₁::Dual, h₂::Dual) = realpart(h₁) == realpart(h₂)
Base.:(==)(h::Dual, x::Number) = realpart(h) == x
Base.:(==)(x::Number, h::Dual) = h == x

Base.isequal(h₁::Dual, h₂::Dual) = isequal(realpart(h₁), realpart(h₂)) && isequal(εpart(h₁), εpart(h₂))
Base.isequal(h::Dual, x::Number) = isequal(realpart(h), x) && isequal(εpart(h), zero(x))
Base.isequal(x::Number, h::Dual) = isequal(h, x)

Base.isless(h₁::Dual{T}, h₂::Dual{T}) where {T <: Real} = realpart(h₁) < realpart(h₂)
Base.isless(h₁::Real, h₂::Dual{<:Real}) = h₁ < realpart(h₂)
Base.isless(h₁::Dual{<:Real}, h₂::Real) = realpart(h₁) < h₂

function Base.hash(h::Dual) # Maybe this works
    x = hash(realpart(h))
    if isequal(h, realpart(h))
        return x
    else
        y = hash(εpart(h))
        return hash(x, hash(y))
    end
end

Base.float(h::Union{Dual{T}, Dual{Complex{T}}}) where {T<:AbstractFloat} = h
Base.complex(h::Dual{<:Complex}) = h
Base.complex(::Type{Dual{T}}) where {T} = Dual{complex(T)}

Base.floor(h::Dual) = floor(realpart(h))
Base.ceil(h::Dual)  = ceil(realpart(h))
Base.trunc(h::Dual) = trunc(realpart(h))
Base.round(h::Dual) = round(realpart(h))
Base.floor(::Type{T}, h::Dual) where {T<:Real} = floor(T, realpart(h))
Base.ceil( ::Type{T}, h::Dual) where {T<:Real} =  ceil(T, realpart(h))
Base.trunc(::Type{T}, h::Dual) where {T<:Real} = trunc(T, realpart(h))
Base.round(::Type{T}, h::Dual) where {T<:Real} = round(T, realpart(h))

for op in (:real, :imag, :conj, :float, :complex)
    @eval Base.$op(z::Dual) = Dual($op(realpart(z)), $op(εpart(z)))
end

Base.abs(z::Dual) = sqrt(abs2(z))
Base.abs2(z::Dual) = real(conj(z) * z)

Base.real(z::Dual{<:Real}) = z
Base.abs(z::Dual{<:Real}) = z ≥ 0 ? z : -z

Base.angle(z::Dual{<:Real}) = z ≥ 0 ? zero(z) : one(z) * π
function Base.angle(z::Dual{Complex{T}}) where T <: Real
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

Base.flipsign(x::Dual, y::Dual) = y == 0 ? flipsign(x, εpart(y)) : flipsign(x, realpart(y))
Base.flipsign(x, y::Dual) = y == 0 ? flipsign(x, εpart(y)) : flipsign(x, realpart(y))
Base.flipsign(x::Dual, y) = dual(flipsign(realpart(x), y), flipsign(εpart(x), y))

Base.:+(h₁::Dual, h₂::Dual) = Dual(realpart(h₁) + realpart(h₂), εpart(h₁) + εpart(h₂))
Base.:+(n::Number, h::Dual) = Dual(n + realpart(h), εpart(h))
Base.:+(h::Dual, n::Number) = n + h

Base.:-(h::Dual) = Dual(-realpart(h), -εpart(h))
Base.:-(h₁::Dual, h₂::Dual) = Dual(realpart(h₁) - realpart(h₂), εpart(h₁) - εpart(h₂))
Base.:-(n::Number, h::Dual) = Dual(n - realpart(h), -εpart(h))
Base.:-(h::Dual, n::Number) = Dual(realpart(h) - n, εpart(h))

# avoid ambiguous definition
Base.:*(x::Bool, h::Dual) = ifelse(x, h, ifelse(signbit(real(realpart(h))) == 0, zero(h), -zero(h)))
Base.:*(h::Dual, x::Bool) = x * h

function Base.:*(h₁::Dual, h₂::Dual)
    x, y = realpart(h₁), εpart(h₁)
    a, b = realpart(h₂), εpart(h₂)
    return Dual(a*x, muladd(a, y, b*x))
end
Base.:*(n::Number, h::Dual) = Dual(n*realpart(h), n*εpart(h))
Base.:*(h::Dual, n::Number) = n * h

Base.one(h::Dual) = Dual(one(realpart(h)))

@inline Base.literal_pow(::typeof(^), x::Dual, ::Val{0}) = one(typeof(x))
@inline Base.literal_pow(::typeof(^), x::Dual, ::Val{1}) = x
@inline Base.literal_pow(::typeof(^), x::Dual, ::Val{2}) = x * x
@inline Base.literal_pow(::typeof(^), x::Dual, ::Val{3}) = x * x * x

function Base.:/(h₁::Dual, h₂::Dual)
    x, y = realpart(h₁), εpart(h₁)
    a, b = realpart(h₂), εpart(h₂)
    return Dual(x / a, y / a - b * x / a ^ 2)
end
function Base.:/(n::Number, h::Dual)
    x, y = realpart(h), ε₁part(h)
    return Dual(n / x, -n * y / x ^ 2)
end
Base.:/(h::Dual, n::Number) = Dual(realpart(h) / n, εpart(h) / n)

Base.mod(h::Dual, n::Number) = Dual(mod(realpart(h), n), εpart(h))

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

Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, i::Dual) = checkindex(Bool, inds, realpart(h))