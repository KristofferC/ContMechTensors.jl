using ForwardDiff
using ContMechTensors
import ForwardDiff: Dual, partials

@inline function _extract{N, T}(v::Dual{N, T})
    Vec{N, T}(partials(v).values)
end

@inline function _load{T}(v::Vec{1, T})
    @inbounds v_dual = Vec{1}((Dual(v[1], one(T)),))
    v_dual
end

@inline function _load{T}(v::Vec{2, T})
    o = one(T)
    z = zero(T)
    @inbounds v_dual = Vec{2}((Dual(v[1], o, z),
                               Dual(v[2], z, o)))
    return v_dual
end

@inline function _load{T}(v::Vec{3, T})
    o = one(T)
    z = zero(T)
    @inbounds v_dual = Vec{3}((Dual(v[1], o, z, z),
                               Dual(v[2], z, o, z),
                               Dual(v[3], z, z, o)))
    return v_dual
end

function gradient{F}(f::F, v::Vec)
    v_dual = _load(v)
    res = f(v_dual)
    return _extract(res)
end

#=
function _extract{D <: Dual}(v::Vec{1, D})
    p1 = partials(v[1])
    Tensor{2, 1}((p1[1],))
end

function _extract{D <: Dual}(v::Vec{2, D})
    p1 = partials(v[1])
    p2 = partials(v[2])
    return Tensor{2, 2}((p1[1], p2[1], p1[2], p2[2]))
end
=#
