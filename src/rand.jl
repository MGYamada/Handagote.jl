Base.rand(r::AbstractRNG, ::Random.SamplerType{Dual{T}}) where T = Dual(rand(r, T), rand(r, T))
Base.rand(r::AbstractRNG, ::Random.SamplerType{HyperDual{T}}) where T = HyperDual(rand(r, T), rand(r, T), rand(r, T), rand(r, T))
Base.rand(r::AbstractRNG, ::Type{Dual{T}}, dims::Tuple{Vararg{Int, N}}) where {N, T} = DualArray(rand(r, T, dims), rand(r, T, dims))
Base.rand(r::AbstractRNG, ::Type{HyperDual{T}}, dims::Tuple{Vararg{Int, N}}) where {N, T} = HyperDualArray(rand(r, T, dims), rand(r, T, dims), rand(r, T, dims), rand(r, T, dims))