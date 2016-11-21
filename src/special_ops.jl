# specialized methods
"""
Computes a special dot product between two vectors and a symmetric fourth order tensor
such that ``a_k C_{ikjl} b_l``.

```julia
dotdot(::Vec, ::SymmetricFourthOrderTensor, ::Vec)
```
"""
@generated function dotdot{dim}(v1::Vec{dim}, S::SymmetricTensor{4, dim}, v2::Vec{dim})
    idx(i,j,k,l) = compute_index(SymmetricTensor{4, dim}, i, j, k, l)
    exps = Expr(:tuple)
    for j in 1:dim, i in 1:dim
        exps_ele = Expr[]
        for l in 1:dim, k in 1:dim
            push!(exps_ele, :(get_data(v1)[$k] * get_data(S)[$(idx(i,k,j,l))] * get_data(v2)[$l]))
        end
        push!(exps.args, reduce((ex1, ex2) -> :(+($ex1, $ex2)), exps_ele))
    end
    return quote
        $(Expr(:meta, :inline))
        @inbounds r = $exps
        Tensor{2, dim}(r)
    end
end

const VOIGT_ORDER = ([1;],
                     [1 3;
                      4  2],
                     [1   6  5;
                      9   2  4;
                      8   7  3])

###########
# tovoigt #
###########
function tovoigt{dim, T}(t::Tensor{2, dim, T})
    v = zeros(T, dim^2)
    for i in 1:dim, j in 1:dim
        v[VOIGT_ORDER[dim][i, j]] = t[i,j]
    end
    return v
end

function tovoigt{dim, T}(t::Tensor{4, dim, T})
    v = zeros(T, dim^2, dim^2)
    for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim
        v[VOIGT_ORDER[dim][i, j], VOIGT_ORDER[dim][k, l]] = t[i,j,k,l]
    end
    return v
end

function tovoigt{dim, T}(t::SymmetricTensor{2, dim, T})
    v = zeros(T, n_components(SymmetricTensor{2, dim}))
    for i in 1:dim, j in i:dim
        v[VOIGT_ORDER[dim][i, j]] = t[i,j]
    end
    return v
end

function tovoigt{dim, T}(t::SymmetricTensor{4, dim, T})
    nc = n_components(SymmetricTensor{2, dim})
    v = zeros(T, nc, nc)
    for i in 1:dim, j in i:dim, k in 1:dim, l in k:dim
        v[VOIGT_ORDER[dim][i, j], VOIGT_ORDER[dim][k, l]] = t[i,j,k,l]
    end
    return v
end
