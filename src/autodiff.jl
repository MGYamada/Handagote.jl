# ChainRules

function ChainRulesCore.rrule(::Type{DualArray}, x, y)
    function dualarray_pullback(ΔΩ)
        (NoTangent(), ΔΩ.value, ΔΩ.epsilon)
    end
    DualArray(x, y), dualarray_pullback
end

function ChainRulesCore.rrule(::typeof(realpart), x::Dual)
    function realpart_pullback(ΔΩ)
        (NoTangent(), (; value = ΔΩ, epsilon = ZeroTangent()))
    end
    realpart(x), realpart_pullback
end

function ChainRulesCore.rrule(::typeof(realpart), x::DualArray)
    function realpart_pullback(ΔΩ)
        (NoTangent(), (; value = ΔΩ, epsilon = ZeroTangent()))
    end
    realpart(x), realpart_pullback
end

function ChainRulesCore.rrule(::typeof(εpart), x::Dual)
    function εpart_pullback(ΔΩ)
        (NoTangent(), (; value = ZeroTangent(), epsilon = ΔΩ))
    end
    εpart(x), εpart_pullback
end

function ChainRulesCore.rrule(::typeof(εpart), x::DualArray)
    function εpart_pullback(ΔΩ)
        (NoTangent(), (; value = ZeroTangent(), epsilon = ΔΩ))
    end
    εpart(x), εpart_pullback
end

function ChainRulesCore.rrule(::Type{HyperDualArray}, x, y, z, w)
    function hyperdualarray_pullback(ΔΩ)
        (NoTangent(), ΔΩ.value, ΔΩ.epsilon1, ΔΩ.epsilon2, ΔΩ.epsilon12)
    end
    HyperDualArray(x, y, z, w), hyperdualarray_pullback
end

function ChainRulesCore.rrule(::typeof(hyperrealpart), x::HyperDual)
    function hyperrealpart_pullback(ΔΩ)
        (NoTangent(), (; value = ΔΩ, epsilon1 = ZeroTangent(),  epsilon2 = ZeroTangent(), epsilon12 = ZeroTangent()))
    end
    hyperrealpart(x), hyperrealpart_pullback
end

function ChainRulesCore.rrule(::typeof(hyperrealpart), x::HyperDualArray)
    function hyperrealpart_pullback(ΔΩ)
        (NoTangent(), (; value = ΔΩ, epsilon1 = ZeroTangent(),  epsilon2 = ZeroTangent(), epsilon12 = ZeroTangent()))
    end
    hyperrealpart(x), hyperrealpart_pullback
end

function ChainRulesCore.rrule(::typeof(ɛ₁part), x::HyperDual)
    function ɛ₁part_pullback(ΔΩ)
        (NoTangent(), (; value = ZeroTangent(), epsilon1 = ΔΩ, epsilon2 = ZeroTangent(), epsilon12 = ZeroTangent()))
    end
    ɛ₁part(x), ɛ₁part_pullback
end

function ChainRulesCore.rrule(::typeof(ɛ₁part), x::HyperDualArray)
    function ɛ₁part_pullback(ΔΩ)
        (NoTangent(), (; value = ZeroTangent(), epsilon1 = ΔΩ, epsilon2 = ZeroTangent(), epsilon12 = ZeroTangent()))
    end
    ɛ₁part(x), ɛ₁part_pullback
end

function ChainRulesCore.rrule(::typeof(ɛ₂part), x::HyperDual)
    function ɛ₂part_pullback(ΔΩ)
        (NoTangent(), (; value = ZeroTangent(), epsilon1 = ZeroTangent(), epsilon2 = ΔΩ, epsilon12 = ZeroTangent()))
    end
    ɛ₂part(x), ɛ₂part_pullback
end

function ChainRulesCore.rrule(::typeof(ɛ₂part), x::HyperDualArray)
    function ɛ₂part_pullback(ΔΩ)
        (NoTangent(), (; value = ZeroTangent(), epsilon1 = ZeroTangent(), epsilon2 = ΔΩ, epsilon12 = ZeroTangent()))
    end
    ɛ₂part(x), ɛ₂part_pullback
end

function ChainRulesCore.rrule(::typeof(ɛ₁ε₂part), x::HyperDual)
    function ɛ₁ε₂part_pullback(ΔΩ)
        (NoTangent(), (; value = ZeroTangent(), epsilon1 = ZeroTangent(), epsilon2 = ZeroTangent(), epsilon12 = ΔΩ))
    end
    ɛ₁ε₂part(x), ɛ₁ε₂part_pullback
end

function ChainRulesCore.rrule(::typeof(ɛ₁ε₂part), x::HyperDualArray)
    function ɛ₁ε₂part_pullback(ΔΩ)
        (NoTangent(), (; value = ZeroTangent(), epsilon1 = ZeroTangent(), epsilon2 = ZeroTangent(), epsilon12 = ΔΩ))
    end
    ɛ₁ε₂part(x), ɛ₁ε₂part_pullback
end

# Zygote

Zygote.@adjoint function LinearAlgebra.tr(x::DualMatrix)
    tr(x), function (Δ)
        ((; value = Diagonal(Fill(Δ.value, (size(x, 1), ))), epsilon = Diagonal(Fill(Δ.epsilon, (size(x, 1), )))),)
    end
end

Zygote.@adjoint function LinearAlgebra.tr(x::HyperDualMatrix)
    tr(x), function (Δ)
        ((; value = Diagonal(Fill(Δ.value, (size(x, 1), ))), epsilon1 = Diagonal(Fill(Δ.epsilon1, (size(x, 1), ))), epsilon2 = Diagonal(Fill(Δ.epsilon2, (size(x, 1), ))), epsilon12 = Diagonal(Fill(Δ.epsilon12, (size(x, 1), )))),)
    end
end