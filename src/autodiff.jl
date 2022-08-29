function ChainRulesCore.rrule(::typeof(tr), x::DualMatrix{T}) where T
    d = size(x, 1)
    function tr_pullback(ΔΩ)
        return (NoTangent(), Tangent{DualMatrix{T}}(; value = Diagonal(fill(ΔΩ.value, d)), epsilon = Diagonal(fill(ΔΩ.epsilon, d))))
    end
    tr(x), tr_pullback
end

function ChainRulesCore.rrule(::typeof(tr), x::HyperDualMatrix{T}) where T
    d = size(x, 1)
    function tr_pullback(ΔΩ)
        return (NoTangent(), Tangent{HyperDualMatrix{T}}(; value = Diagonal(fill(ΔΩ.value, d)), epsilon1 = Diagonal(fill(ΔΩ.epsilon1, d)), epsilon2 = Diagonal(fill(ΔΩ.epsilon2, d)), epsilon12 = Diagonal(fill(ΔΩ.epsilon12, d))))
    end
    tr(x), tr_pullback
end