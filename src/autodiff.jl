Zygote.@adjoint function Dual{T}(x, y) where T
    Dual{T}(x, y), Δ -> (Δ.value, Δ.epsilon)
end

Zygote.@adjoint function DualArray{T, N}(x, y) where {T, N}
    DualArray{T, N}(x, y), Δ -> (Δ.value, Δ.epsilon)
end

Zygote.@adjoint function realpart(x::Dual)
    realpart(x), Δ -> (Dual(Δ, zero(Δ)),)
end

Zygote.@adjoint function realpart(x::DualArray)
    realpart(x), Δ -> (DualArray(Δ, zero.(Δ)),)
end

Zygote.@adjoint function εpart(x::Dual)
    εpart(x), Δ -> (Dual(zero(Δ), Δ),)
end

Zygote.@adjoint function εpart(x::DualArray)
    εpart(x), Δ -> (DualArray(zero.(Δ), Δ),)
end

Zygote.@adjoint function HyperDual{T}(x, y, z, w) where T
    HyperDual{T}(x, y, z, w), Δ -> (Δ.value, Δ.epsilon1, Δ.epsilon2, Δ.epsilon12)
end

Zygote.@adjoint function HyperDualArray{T, N}(x, y, z, w) where {T, N}
    HyperDualArray{T, N}(x, y, z, w), Δ -> (Δ.value, Δ.epsilon1, Δ.epsilon2, Δ.epsilon12)
end

Zygote.@adjoint function hyperrealpart(x::HyperDual)
    hyperrealpart(x), Δ -> (HyperDual(Δ, zero(Δ), zero(Δ), zero(Δ)),)
end

Zygote.@adjoint function hyperrealpart(x::HyperDualArray)
    hyperrealpart(x), Δ -> (HyperDualArray(Δ, zero.(Δ), zero.(Δ), zero.(Δ)),)
end

Zygote.@adjoint function ɛ₁part(x::HyperDual)
    ɛ₁part(x), Δ -> (HyperDual(zero(Δ), Δ, zero(Δ), zero(Δ)),)
end

Zygote.@adjoint function ɛ₁part(x::HyperDualArray)
    ɛ₁part(x), Δ -> (HyperDualArray(zero.(Δ), Δ, zero.(Δ), zero.(Δ)),)
end

Zygote.@adjoint function ɛ₂part(x::HyperDual)
    ɛ₂part(x), Δ -> (HyperDual(zero(Δ), zero(Δ), Δ, zero(Δ)),)
end

Zygote.@adjoint function ɛ₂part(x::HyperDualArray)
    ɛ₂part(x), Δ -> (HyperDualArray(zero.(Δ), zero.(Δ), Δ, zero.(Δ)),)
end

Zygote.@adjoint function ɛ₁ε₂part(x::HyperDual)
    ɛ₁ε₂part(x), Δ -> (HyperDual(zero(Δ), zero(Δ), zero(Δ), Δ),)
end

Zygote.@adjoint function ɛ₁ε₂part(x::HyperDualArray)
    ɛ₁ε₂part(x), Δ -> (HyperDualArray(zero.(Δ), zero.(Δ), zero.(Δ), Δ),)
end

Zygote.@adjoint function LinearAlgebra.tr(x::DualMatrix)
    tr(x), Δ -> (DualArray(Diagonal(Fill(Δ.value, (size(x, 1),))), Diagonal(Fill(Δ.epsilon, (size(x, 1),)))),)
end

Zygote.@adjoint function LinearAlgebra.tr(x::HyperDualMatrix)
    tr(x), Δ -> (HyperDualArray(Diagonal(Fill(Δ.value, (size(x, 1),))), Diagonal(Fill(Δ.epsilon1, (size(x, 1),))), Diagonal(Fill(Δ.epsilon2, (size(x, 1),))), Diagonal(Fill(Δ.epsilon12, (size(x, 1),)))),)
end

Zygote.@adjoint function Base.reshape(xs::DualArray, dims...)
    reshape(xs, dims...), Δ -> (DualArray(reshape(Δ.value, size(xs)), reshape(Δ.epsilon, size(xs))), map(_ -> nothing, dims)...)
end

Zygote.@adjoint function Base.reshape(xs::HyperDualArray, dims...)
    reshape(xs, dims...), Δ -> (HyperDualArray(reshape(Δ.value, size(xs)), reshape(Δ.epsilon1, size(xs)), reshape(Δ.epsilon2, size(xs)), reshape(Δ.epsilon12, size(xs))), map(_ -> nothing, dims)...)
end