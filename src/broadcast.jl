struct DualArrayStyle{N} <: Broadcast.AbstractArrayStyle{N} end
DualArrayStyle(::Val{N}) where N = DualArrayStyle{N}()

Base.BroadcastStyle(::DualArrayStyle{N}, ::Base.Broadcast.DefaultArrayStyle{M}) where {M, N} = DualArrayStyle(Val(max(M, N)))
Base.BroadcastStyle(::DualArrayStyle{N}, ::Base.Broadcast.AbstractArrayStyle{M}) where {M, N} = DualArrayStyle(Val(max(M, N)))
Base.BroadcastStyle(::DualArrayStyle{M}, ::DualArrayStyle{N}) where {M, N} = DualArrayStyle(Val(max(M, N)))
Base.BroadcastStyle(::Type{<:DualArray{T, N}}) where {T, N} = DualArrayStyle{N}()

Base.similar(bc::Broadcast.Broadcasted{DualArrayStyle{N}}, ::Type{T}) where {N, T} = zeros(T, axes(bc))

struct HyperDualArrayStyle{N} <: Broadcast.AbstractArrayStyle{N} end
HyperDualArrayStyle(::Val{N}) where N = HyperDualArrayStyle{N}()

Base.BroadcastStyle(::HyperDualArrayStyle{N}, ::Base.Broadcast.DefaultArrayStyle{M}) where {M, N} = HyperDualArrayStyle(Val(max(M, N)))
Base.BroadcastStyle(::HyperDualArrayStyle{N}, ::Base.Broadcast.AbstractArrayStyle{M}) where {M, N} = HyperDualArrayStyle(Val(max(M, N)))
Base.BroadcastStyle(::HyperDualArrayStyle{M}, ::DualArrayStyle{N}) where {M, N} = HyperDualArrayStyle(Val(max(M, N)))
Base.BroadcastStyle(::HyperDualArrayStyle{M}, ::HyperDualArrayStyle{N}) where {M, N} = HyperDualArrayStyle(Val(max(M, N)))
Base.BroadcastStyle(::Type{<:HyperDualArray{T, N}}) where {T, N} = HyperDualArrayStyle{N}()

Base.similar(bc::Broadcast.Broadcasted{HyperDualArrayStyle{N}}, ::Type{T}) where {N, T} = zeros(T, axes(bc))