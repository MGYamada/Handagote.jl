struct HyperDual{T <: ReComp} <: Number
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

Base.convert(::Type{T}, h::HyperDual) where {T <: Dual} = (ɛ₂part(h) == 0 && ɛ₁ε₂part(h) == 0 ? convert(T, hyperrealpart(h)) : throw(InexactError()))
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

# show

Base.convert(::Type{HyperDual}, h::HyperDual) = h
Base.convert(::Type{HyperDual}, x::Number) = HyperDual(x)

Base.:(==)(h₁::HyperDual, h₂::HyperDual) = hyperrealpart(h₁) == hyperrealpart(h₂)
Base.:(==)(h::HyperDual, x::Number) = hyperrealpart(h) == x
Base.:(==)(x::Number, h::HyperDual) = h == x

Base.isequal(h₁::HyperDual, h₂::HyperDual) = isequal(hyperrealpart(h₁),hyperrealpart(h₂)) && isequal(ε₁part(h₁), ε₁part(h₂)) && isequal(ε₂part(h₁), ε₂part(h₂)) && isequal(ε₁ε₂part(h₁), ε₁ε₂part(h₂))
Base.isequal(h::HyperDual, x::Number) = isequal(hyperrealpart(h), x) && isequal(ε₁part(h), zero(x)) && isequal(ε₂part(h), zero(x)) && isequal(ε₁ε₂part(h), zero(x))
Base.isequal(x::Number, h::HyperDual) = isequal(h, x)

Base.isless(h₁::HyperDual{T}, h₂::HyperDual{T}) where {T<:Real} = hyperrealpart(h₁) < hyperrealpart(h₂)
Base.isless(h₁::Real, h₂::HyperDual{<:Real}) = h₁ < hyperrealpart(h₂)
Base.isless(h₁::HyperDual{<:Real}, h₂::Real) = hyperrealpart(h₁) < h₂

function Base.hash(h::HyperDual) # Not sure this works
    x = hash(hyperrealpart(h))
    if isequal(h, hyperrealpart(h))
        return x
    else
        y = hash(ε₁part(h))
        z = hash(ε₂part(h))
        w = hash(ε₁ε₂part(h))
        return hash(x, hash(y, hash(z, hash(w))))
    end
end

Base.float(h::Union{HyperDual{T}, HyperDual{Complex{T}}}) where {T<:AbstractFloat} = h
Base.complex(h::HyperDual{<:Complex}) = h
Base.complex(::Type{HyperDual{T}}) where {T} = HyperDual{complex(T)}

Base.floor(h::HyperDual) = floor(hyperrealpart(h))
Base.ceil(h::HyperDual)  = ceil(hyperrealpart(h))
Base.trunc(h::HyperDual) = trunc(hyperrealpart(h))
Base.round(h::HyperDual) = round(hyperrealpart(h))
Base.floor(::Type{T}, h::HyperDual) where {T<:Real} = floor(T, hyperrealpart(h))
Base.ceil( ::Type{T}, h::HyperDual) where {T<:Real} =  ceil(T, hyperrealpart(h))
Base.trunc(::Type{T}, h::HyperDual) where {T<:Real} = trunc(T, hyperrealpart(h))
Base.round(::Type{T}, h::HyperDual) where {T<:Real} = round(T, hyperrealpart(h))

for op in (:real, :imag, :conj, :float, :complex)
    @eval Base.$op(z::HyperDual) = HyperDual($op(realpart(z)), $op(εpart(z)))
end

Base.abs(z::HyperDual) = sqrt(abs2(z))
Base.abs2(z::HyperDual) = real(conj(z) * z)

Base.real(z::HyperDual{<:Real}) = z
Base.abs(z::HyperDual{<:Real}) = z ≥ 0 ? z : -z

Base.angle(z::HyperDual{<:Real}) = z ≥ 0 ? zero(z) : one(z) * π
function Base.angle(z::HyperDual{Complex{T}}) where T <: Real
    if z == 0
        HyperDual(zero(T), imag(ɛ₁part(z)) == 0 ? zero(T) : convert(T, Inf), imag(ε₂part(z)) == 0 ? zero(T) : convert(T, Inf), imag(ɛ₁ε₂part(z)) == 0 ? zero(T) : convert(T, Inf))
    else
        real(log(sign(z)) / im)
    end
end

# Base.flipsign is undefined

Base.:+(h₁::HyperDual, h₂::HyperDual) = HyperDual(hyperrealpart(h₁) + hyperrealpart(h₂), ε₁part(h₁) + ε₁part(h₂), ε₂part(h₁) + ε₂part(h₂), ε₁ε₂part(h₁) + ε₁ε₂part(h₂))
Base.:+(n::Number, h::HyperDual) = HyperDual(n + hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h))
Base.:+(h::HyperDual, n::Number) = n + h

Base.:-(h::HyperDual) = HyperDual(-hyperrealpart(h), -ε₁part(h), -ε₂part(h), -ε₁ε₂part(h))
Base.:-(h₁::HyperDual, h₂::HyperDual) = HyperDual(hyperrealpart(h₁) - hyperrealpart(h₂), ε₁part(h₁) - ε₁part(h₂), ε₂part(h₁) - ε₂part(h₂), ε₁ε₂part(h₁) - ε₁ε₂part(h₂))
Base.:-(n::Number, h::HyperDual) = HyperDual(n - hyperrealpart(h), -ε₁part(h), -ε₂part(h), -ε₁ε₂part(h))
Base.:-(h::HyperDual, n::Number) = HyperDual(hyperrealpart(h) - n, ε₁part(h), ε₂part(h), ε₁ε₂part(h))

# avoid ambiguous definition with Bool*Number
Base.:*(x::Bool, h::HyperDual) = ifelse(x, h, ifelse(signbit(real(hyperrealpart(h)))==0, zero(h), -zero(h)))
Base.:*(h::HyperDual, x::Bool) = x * h

function Base.:*(h₁::HyperDual, h₂::HyperDual)
    x, y, z, w = hyperrealpart(h₁), ε₁part(h₁), ε₂part(h₁), ε₁ε₂part(h₁)
    a, b, c, d = hyperrealpart(h₂), ε₁part(h₂), ε₂part(h₂), ε₁ε₂part(h₂)
    return HyperDual(a*x, muladd(a, y, b*x), muladd(a, z, c*x), muladd(a, w, muladd(d, x, muladd(c, y, b*z))))
end
Base.:*(n::Number, h::HyperDual) = HyperDual(n*hyperrealpart(h), n*ε₁part(h), n*ε₂part(h), n*ε₁ε₂part(h))
Base.:*(h::HyperDual, n::Number) = n * h

Base.one(h::HyperDual) = HyperDual(one(realpart(h)))

@inline Base.literal_pow(::typeof(^), x::HyperDual, ::Val{0}) = one(typeof(x))
@inline Base.literal_pow(::typeof(^), x::HyperDual, ::Val{1}) = x
@inline Base.literal_pow(::typeof(^), x::HyperDual, ::Val{2}) = x*x
@inline Base.literal_pow(::typeof(^), x::HyperDual, ::Val{3}) = x*x*x

function Base.:/(h₁::HyperDual, h₂::HyperDual)
    x, y, z, w = hyperrealpart(h₁), ε₁part(h₁), ε₂part(h₁), ε₁ε₂part(h₁)
    a, b, c, d = hyperrealpart(h₂), ε₁part(h₂), ε₂part(h₂), ε₁ε₂part(h₂)
    return HyperDual(x/a, y/a - b*x/a^2, z/a - c*x/a^2, w/a + 2*b*c*x/a^3 - d*x/a^2 - c*y/a^2 - b*z/a^2)
end
function Base.:/(n::Number, h::HyperDual)
    x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    return HyperDual(n/x, -n*y/x^2, -n*z/x^2, -n*(w/x-y*z/x^2-z*y/x^2)/x)
end
Base.:/(h::HyperDual, n::Number) = HyperDual(hyperrealpart(h)/n, ε₁part(h)/n, ε₂part(h)/n, ε₁ε₂part(h)/n)

Base.mod(h::HyperDual, n::Number) = Hyper(mod(hyperrealpart(h), n), ε₁part(h), ε₂part(h), ε₁ε₂part(h))

# pow

for (fsym, dfexp, d²fexp) in symbolic_derivative_list
    mod = isdefined(SpecialFunctions, fsym) ? SpecialFunctions :
          isdefined(Base, fsym)             ? Base             :
          isdefined(Base.Math, fsym)        ? Base.Math        :
          nothing
    if mod !== nothing && fsym !== :sin && fsym !== :cos # (we define out own sin and cos)
        expr = :(HyperDual($(fsym)(x), y*$dfexp, z*$dfexp, w*$dfexp + y*z*$d²fexp))
        cse_expr = CommonSubexpressions.cse(expr, warn=false)

        @eval function $mod.$(fsym)(h::HyperDual)
            x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
            $cse_expr
        end
    end
    # extend corresponding NaNMath methods
    if fsym in (:sin, :cos, :tan, :asin, :acos, :acosh, :atanh, :log, :log2, :log10, :log1p)
        fsym = Expr(:.,:NaNMath,Base.Meta.quot(fsym))
        @eval function $(fsym)(h::HyperDual)
            x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
            HyperDual($(fsym)(x), y*$(to_nanmath(dfexp)), z*$(to_nanmath(dfexp)), w*$(to_nanmath(dfexp)) + y*z*$(to_nanmath(d²fexp)))
        end
    end
end