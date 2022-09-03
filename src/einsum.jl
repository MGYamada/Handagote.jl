macro dualein_str(s::AbstractString)
    dualein(s)
end

struct DualStaticEinCode{LT}
    code::StaticEinCode{LT}
end

code(c::DualStaticEinCode) = c.code

function dualein(s::AbstractString)
    DualStaticEinCode(ein(s))
end

function (c::DualStaticEinCode{LT})(xs...; kwargs...) where LT
    apply_linear(code(c), xs...; kwargs...)
end

macro hyperdualein_str(s::AbstractString)
    hyperdualein(s)
end

struct HyperDualStaticEinCode{LT}
    code::DualStaticEinCode{LT}
end

code(c::HyperDualStaticEinCode) = c.code

function hyperdualein(s::AbstractString)
    HyperDualStaticEinCode(dualein(s))
end

function (c::HyperDualStaticEinCode{LT})(xs...; kwargs...) where LT
    apply_linear_hyper(code(c), xs...; kwargs...)
end