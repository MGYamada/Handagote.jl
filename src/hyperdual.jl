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

function hyperpart_show(io::IO, yzw::T, compact::Bool, str::String) where T<:Real
    if signbit(yzw)
        yzw = -yzw
        print(io, compact ? "-" : " - ")
    else
        print(io, compact ? "+" : " + ")
    end
    compact ? show(IOContext(io, :compact=>true), yzw) : show(io, yzw)
    printtimes(io, yzw)
    print(io, str)
end

function hyper_show(io::IO, h::HyperDual{T}, compact::Bool) where T<:Real
    x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    compact ? show(IOContext(io, :compact=>true), x) : show(io, x)
    hyperpart_show(io, y, compact, "ε₁")
    hyperpart_show(io, z, compact, "ε₂")
    hyperpart_show(io, w, compact, "ε₁ε₂")
end

function hyperpart_show(io::IO, yzw::T, compact::Bool, str::String) where T<:Complex
    yzwr, yzwi = reim(yzw)
    if signbit(yzwr)
        yzwr = -yzwr
        print(io, " - ")
    else
        print(io, " + ")
    end
    if compact
        if signbit(yzwi)
            yzwi = -yzwi
            show(IOContext(io, :compact=>true), yzwr)
            printtimes(io, yzwr)
            print(io, str, "-")
            show(IOContext(io, :compact=>true), yzwi)
        else
            show(IOContext(io, :compact=>true), yzwr)
            print(io, str, "+")
            show(IOContext(io, :compact=>true), yzwi)
        end
    else
        if signbit(yzwi)
            yzwi = -yzwi
            show(io, yzwr)
            printtimes(io, yzwr)
            print(io, str, " - ")
            show(io, yzwi)
        else
            show(io, yzwr)
            print(io, str, " + ")
            show(io, yzwi)
        end
    end
    printtimes(io, yzwi)
    print(io, "im", str)
end

function hyper_show(io::IO, h::HyperDual{T}, compact::Bool) where T<:Complex
    x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    compact ? show(IOContext(io, :compact=>true), x) : show(io, x)
    hyperpart_show(io, y, compact, "ε₁")
    hyperpart_show(io, z, compact, "ε₂")
    hyperpart_show(io, w, compact, "ε₁ε₂")
end

function hyper_show(io::IO, h::HyperDual{T}, compact::Bool) where T<:Bool
    x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    if !x && y && !z && !w
        print(io, "ɛ₁")
    elseif !x && !y && z && !w
        print(io, "ɛ₂")
    elseif !x && !y && !z && w
        print(io, "ɛ₁ε₂")
    else
        print(io, "HyperDual{",T,"}(", x, ",", y, ",", z, ",", w, ")")
    end
end

function hyper_show(io::IO, h::HyperDual{Complex{T}}, compact::Bool) where T<:Bool
    x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    xr, xi = reim(x)
    yr, yi = reim(y)
    zr, zi = reim(z)
    wr, wi = reim(w)
    if !xr * xi * !yr * !yi * !zr * !zi * !wr * !wi
        print(io, "im")
    elseif !xr * !xi * yr * !yi * !zr * !zi * !wr * !wi
        print(io, "ɛ₁")
    elseif !xr * !xi * !yr * yi * !zr * !zi * !wr * !wi
        print(io, "imɛ₁")
    elseif !xr * !xi * !yr * !yi * zr * !zi * !wr * !wi
        print(io, "ɛ₂")
    elseif !xr * !xi * !yr * !yi * !zr * zi * !wr * !wi
        print(io, "imɛ₂")
    elseif !xr * !xi * !yr * !yi * !zr * !zi * wr * !wi
        print(io, "ε₁ɛ₂")
    elseif !xr * !xi * !yr * !yi * !zr * !zi * !wr * wi
        print(io, "imε₁ɛ₂")
    else
        print(io, "HyperDual{",T,"}(", x, ",", y, ",", z, ",", w, ")")
    end
end

function printtimes(io::IO, x::Real)
    if !(isa(x,Integer) || isa(x,Rational) ||
         isa(x,AbstractFloat) && isfinite(x))
        print(io, "*")
    end
end

Base.show(io::IO, h::HyperDual) = hyper_show(io, h, get(IOContext(io), :compact, false))

function Base.read(s::IO, ::Type{HyperDual{T}}) where T<:ReComp
    x = read(s, T)
    y = read(s, T)
    z = read(s, T)
    w = read(s, T)
    HyperDual{T}(x, y, z, w)
end
function Base.write(s::IO, h::HyperDual)
    write(s, hyperrealpart(h))
    write(s, ε₁part(h))
    write(s, ε₂part(h))
    write(s, ε₁ε₂part(h))
end

Base.convert(::Type{HyperDual}, h::HyperDual) = h
Base.convert(::Type{HyperDual}, x::Number) = HyperDual(x)

Base.:(==)(h₁::HyperDual, h₂::HyperDual) = hyperrealpart(h₁) == hyperrealpart(h₂)
Base.:(==)(h::HyperDual, x::Number) = hyperrealpart(h) == x
Base.:(==)(x::Number, h::HyperDual) = h == x

Base.isequal(h₁::HyperDual, h₂::HyperDual) = isequal(hyperrealpart(h₁), hyperrealpart(h₂)) && isequal(ε₁part(h₁), ε₁part(h₂)) && isequal(ε₂part(h₁), ε₂part(h₂)) && isequal(ε₁ε₂part(h₁), ε₁ε₂part(h₂))
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

Base.mod(h::HyperDual, n::Number) = HyperDual(mod(hyperrealpart(h), n), ε₁part(h), ε₂part(h), ε₁ε₂part(h))

# Power functions written using sage to see Taylor expansions
#   (x+y*ε₁+z*ε₂+w*ε₁*ε₂)^(a+b*ε₁+c*ε₂+d*ε₁*ε₂)
# around 0 for y, z, w, b, c, and d
function Base.:^(h₁::HyperDual, h₂::HyperDual)
    x, y, z, w = hyperrealpart(h₁), ε₁part(h₁), ε₂part(h₁), ε₁ε₂part(h₁)
    a, b, c, d = hyperrealpart(h₂), ε₁part(h₂), ε₂part(h₂), ε₁ε₂part(h₂)
    return HyperDual(x^a,
        a*x^(a - 1)*y + b*x^a*log(x),
        a*x^(a - 1)*z + c*x^a*log(x),
        a^2*x^(a - 2)*y*z + a*c*x^(a - 1)*y*log(x) + a*b*x^(a - 1)*z*log(x) + b*c*x^a*log(x)^2 - a*x^(a - 2)*y*z + a*w*x^(a - 1) + c*x^(a - 1)*y + b*x^(a - 1)*z + d*x^a*log(x))
end

function NaNMath.pow(h₁::HyperDual, h₂::HyperDual)
    x, y, z, w = hyperrealpart(h₁), ε₁part(h₁), ε₂part(h₁), ε₁ε₂part(h₁)
    a, b, c, d = hyperrealpart(h₂), ε₁part(h₂), ε₂part(h₂), ε₁ε₂part(h₂)
    return HyperDual(NaNMath.pow(x,a),
        a*NaNMath.pow(x,a - 1)*y + b*NaNMath.pow(x,a)*log(x),
        a*NaNMath.pow(x,a - 1)*z + c*NaNMath.pow(x,a)*log(x),
        a^2*NaNMath.pow(x,a - 2)*y*z + a*c*NaNMath.pow(x,a - 1)*y*log(x) + a*b*NaNMath.pow(x,a - 1)*z*log(x) + b*c*NaNMath.pow(x,a)*log(x)^2 - a*NaNMath.pow(x,a - 2)*y*z + a*w*NaNMath.pow(x,a - 1) + c*NaNMath.pow(x,a - 1)*y + b*NaNMath.pow(x,a - 1)*z + d*NaNMath.pow(x,a)*log(x))
end

function Base.:^(h::HyperDual, a::Integer)
    x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    return HyperDual(x^a,
        a*x^(a - 1)*y,
        a*x^(a - 1)*z,
        a^2*x^(a - 2)*y*z - a*x^(a - 2)*y*z + a*w*x^(a - 1))
end

function Base.:^(h::HyperDual, a::Rational)
    x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    return HyperDual(x^a,
        a*x^(a - 1)*y,
        a*x^(a - 1)*z,
        a^2*x^(a - 2)*y*z - a*x^(a - 2)*y*z + a*w*x^(a - 1))
end

function Base.:^(h::HyperDual, a::Number)
    x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    return HyperDual(x^a,
        a*x^(a - 1)*y,
        a*x^(a - 1)*z,
        a^2*x^(a - 2)*y*z - a*x^(a - 2)*y*z + a*w*x^(a - 1))
end

# Below definition is necesssary to resolve a conflict with the
# definition in MathConstants.jl
function Base.:^(x::Irrational{:ℯ}, h::HyperDual)
    a, b, c, d = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    return HyperDual(x^a,
        b*x^a*log(x),
        c*x^a*log(x),
        b*c*x^a*log(x)^2 + d*x^a*log(x))
end
function Base.:^(x::Number, h::HyperDual)
    a, b, c, d = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    return HyperDual(x^a,
        b*x^a*log(x),
        c*x^a*log(x),
        b*c*x^a*log(x)^2 + d*x^a*log(x))
end

function NaNMath.pow(h::HyperDual, a::Number)
    x, y, z, w = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    return HyperDual(NaNMath.pow(x,a),
        a*NaNMath.pow(x,a - 1)*y,
        a*NaNMath.pow(x,a - 1)*z,
        a^2*NaNMath.pow(x,a - 2)*y*z - a*NaNMath.pow(x,a - 2)*y*z + a*w*NaNMath.pow(x,a - 1))
end
function NaNMath.pow(x::Number, h::HyperDual)
    a, b, c, d = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    return HyperDual(NaNMath.pow(x,a),
        b*NaNMath.pow(x,a)*log(x),
        c*NaNMath.pow(x,a)*log(x),
        b*c*NaNMath.pow(x,a)*log(x)^2 + d*NaNMath.pow(x,a)*log(x))
end

# function to_nanmath(x::Expr)
#     if x.head == :call
#         funsym = Expr(:.,:NaNMath,Base.Meta.quot(x.args[1]))
#         return Expr(:call,funsym,[to_nanmath(z) for z in x.args[2:end]]...)
#     else
#         return Expr(:call,[to_nanmath(z) for z in x.args]...)
#     end
# end
# to_nanmath(x) = x

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

# Can use sincos for cos and sin
function Base.cos(h::HyperDual)
    a, b, c, d = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    si, co = sincos(a)
    return HyperDual(co, -si*b, -si*c, -si*d - co*b*c)
end

function Base.sin(h::HyperDual)
    a, b, c, d = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    si, co = sincos(a)
    return HyperDual(si, co*b, co*c, co*d - si*b*c)
end

# only need to compute exp/cis once (removed exp from derivatives_list)
function Base.exp(h::HyperDual)
    a, b, c, d = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    return exp(a) * HyperDual(one(a), b, c, d + b*c)
end

function Base.cis(h::HyperDual)
    a, b, c, d = hyperrealpart(h), ε₁part(h), ε₂part(h), ε₁ε₂part(h)
    return cis(a) * HyperDual(one(a), im*b, im*c, im*d - b*c)
end

# TODO: should be generated in Calculus, sinpi and cospi (erased here)

Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, i::HyperDual) = checkindex(Bool, inds, hyperrealpart(h))