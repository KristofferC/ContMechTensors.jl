#################################################
# Sepcialized Second Order Symmetric Operations #
#################################################

Base.transpose(S::SymmetricTensors) = S
Base.issym(S::SymmetricTensors) = true

######################
# Double contraction #
######################
@gen_code function dcontract{dim, T1, T2}(S1::SymmetricTensor{2, dim, T1}, S2::SymmetricTensor{2, dim, T2})
    Tv = typeof(zero(T1) * zero(T2))
    @code :(s = zero($Tv);
            data1 = get_data(S1);
            data2 = get_data(S2))
     for k in 1:n_independent_components(dim, true)
        if is_diagonal_index(dim, k)
            @code :(@inbounds s += data1[$k] * data2[$k])
        else
            @code :(@inbounds s += 2 * data1[$k] * data2[$k])
        end
    end
    @code :(return s)
end

@generated function dcontract{dim, T1, T2}(S1::SymmetricTensor{2, dim, T1}, S2::SymmetricTensor{4, dim, T2})
    idx4(i,j,k,l) = compute_index(SymmetricTensor{4, dim}, i, j, k, l)
    idx2(k,l) = compute_index(SymmetricTensor{2, dim}, k, l)
    exps = Expr(:tuple)
    for i in 1:dim, j in 1:i
        exps_ele = Expr(:call)
        push!(exps_ele.args, :+)
        for k in 1:dim, l in 1:k
            if k == l
                push!(exps_ele.args, :(data4[$(idx4(k,l,i,j))] * data2[$(idx2(k,l))]))
            else
                push!(exps_ele.args, :( 2 * data4[$(idx4(k,l,i,j))] * data2[$(idx2(k,l))]))
            end
        end
        push!(exps.args, exps_ele)
    end
    quote
         data2 = S1.data
         data4 = S2.data
         @inbounds r = $exps
         SymmetricTensor{2, dim}(r)
    end
end

@generated function dcontract{dim, T1, T2}(S1::SymmetricTensor{4, dim, T1}, S2::SymmetricTensor{2, dim, T2})
    idx4(i,j,k,l) = compute_index(SymmetricTensor{4, dim}, i, j, k, l)
    idx2(k,l) = compute_index(SymmetricTensor{2, dim}, k, l)
    exps = Expr(:tuple)
    for i in 1:dim, j in 1:i
        exps_ele = Expr(:call)
        push!(exps_ele.args, :+)
        for k in 1:dim, l in 1:k
            if k == l
                push!(exps_ele.args, :(data4[$(idx4(i,j,k,l))] * data2[$(idx2(k,l))]))
            else
                push!(exps_ele.args, :( 2 * data4[$(idx4(i,j,k,l))] * data2[$(idx2(k,l))]))
            end
        end
        push!(exps.args, exps_ele)
    end
    quote
         data2 = S2.data
         data4 = S1.data
         @inbounds r = $exps
         SymmetricTensor{2, dim}(r)
    end
end

########
# norm #
########
@gen_code function Base.norm{dim, T}(S::SymmetricTensor{4, dim, T})
    idx(i,j,k,l) = compute_index(SymmetricTensor{4, dim}, i, j, k, l)
    @code :(data = get_data(S))
    @code :(s = zero(T))
    for k in 1:dim, l in 1:k, i in 1:dim, j in 1:i
        @code :(@inbounds v = data[$(idx(i,j,k,l))])
        if i == j && k == l
             @code :(s += v*v)
        elseif i == j || k == l
             @code :(s += 2*v*v)
        else
             @code :(s += 4*v*v)
        end
    end
    @code :(return sqrt(s))
end


################
# Open product #
################
function otimes{dim, T1, T2}(S1::SymmetricTensor{2, dim, T1}, S2::SymmetricTensor{2, dim, T2})
    SymmetricTensor{4, dim}(A_otimes_B(S1.data, S2.data))
end
