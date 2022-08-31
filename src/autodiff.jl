function ChainRulesCore.rrule(::Type{HyperDualArray}, x, y, z, w)
    function hyperdualarray_pullback(ΔΩ)
        (NoTangent(), hyperrealpart(ΔΩ), ɛ₁part(ΔΩ), ɛ₂part(ΔΩ), ɛ₁ε₂part(ΔΩ))
    end
    HyperDualArray(x, y, z, w), hyperdualarray_pullback
end

function ChainRulesCore.rrule(::typeof(hyperrealpart), x::HyperDual)
    function hyperrealpart_pullback(ΔΩ)
        (NoTangent(), HyperDual(ΔΩ))
    end
    hyperrealpart(x), hyperrealpart_pullback
end

function ChainRulesCore.rrule(::typeof(hyperrealpart), x::HyperDualArray)
    function hyperrealpart_pullback(ΔΩ)
        (NoTangent(), HyperDualArray(ΔΩ))
    end
    hyperrealpart(x), hyperrealpart_pullback
end

Zygote.@adjoint function LinearAlgebra.tr(x::HyperDualMatrix)
    tr(x), function (Δ)
        (HyperDualArray(Diagonal(Fill(Δ.value, (size(x, 1), ))), Diagonal(Fill(Δ.epsilon1, (size(x, 1), ))), Diagonal(Fill(Δ.epsilon2, (size(x, 1), ))), Diagonal(Fill(Δ.epsilon12, (size(x, 1), )))),)
    end
end