__precompile__()

module ContMechTensors

import Base.@pure

using StaticArrays

immutable InternalError <: Exception end

export AbstractTensor, SymmetricTensor, Tensor, Vec, FourthOrderTensor, SecondOrderTensor

export otimes, ⊗, ⊡, dcontract, dev, vol, symmetric, skew, minorsymmetric, majorsymmetric
export minortranspose, majortranspose, isminorsymmetric, ismajorsymmetric
export tdot, dotdot

@deprecate extract_components(tensor) Array(tensor)

#########
# Types #
#########

abstract AbstractTensor{order, dim, T <: Real} <: AbstractArray{T, order}

immutable SymmetricTensor{order, dim, T <: Real, M} <: AbstractTensor{order, dim, T}
   data::SVector{M, T}
end
# Tuple constructor
(::Type{SymmetricTensor{order, dim}}){order, dim, T <: Real, M}(t::NTuple{M, T}) = SymmetricTensor{order, dim}(SVector{M,T}(t))

# SVector constructors
(::Type{SymmetricTensor{order, dim}}){order, dim, T <: Real, M}(t::SVector{M, T}) = throw(ArgumentError("error")) # Fallback

(::Type{SymmetricTensor{2, 1}}){T <: Real}(t::SVector{1, T}) = SymmetricTensor{2, 1, T, 1}(t)
(::Type{SymmetricTensor{2, 2}}){T <: Real}(t::SVector{3, T}) = SymmetricTensor{2, 2, T, 3}(t)
(::Type{SymmetricTensor{2, 3}}){T <: Real}(t::SVector{6, T}) = SymmetricTensor{2, 3, T, 6}(t)

(::Type{SymmetricTensor{4, 1}}){T <: Real}(t::SVector{1, T}) = SymmetricTensor{4, 1, T, 1}(t)
(::Type{SymmetricTensor{4, 2}}){T <: Real}(t::SVector{9, T}) = SymmetricTensor{4, 2, T, 9}(t)
(::Type{SymmetricTensor{4, 3}}){T <: Real}(t::SVector{36, T}) = SymmetricTensor{4, 3, T, 36}(t)

immutable Tensor{order, dim, T <: Real, M} <: AbstractTensor{order, dim, T}
   data::SVector{M, T}
end

# Tuple constructor
(::Type{Tensor{order, dim}}){order, dim, T <: Real, M}(t::NTuple{M, T}) = Tensor{order, dim}(SVector{M,T}(t))

# SVector constructors
(::Type{Tensor{order, dim}}){order, dim, T <: Real, M}(t::SVector{M, T}) = throw(ArgumentError("error")) # Fallback

(::Type{Tensor{1, 1}}){T <: Real}(t::SVector{1, T}) = Tensor{1, 1, T, 1}(t)
(::Type{Tensor{1, 2}}){T <: Real}(t::SVector{2, T}) = Tensor{1, 2, T, 2}(t)
(::Type{Tensor{1, 3}}){T <: Real}(t::SVector{3, T}) = Tensor{1, 3, T, 3}(t)

(::Type{Tensor{2, 1}}){T <: Real}(t::SVector{1, T}) = Tensor{2, 1, T, 1}(t)
(::Type{Tensor{2, 2}}){T <: Real}(t::SVector{4, T}) = Tensor{2, 2, T, 4}(t)
(::Type{Tensor{2, 3}}){T <: Real}(t::SVector{9, T}) = Tensor{2, 3, T, 9}(t)

(::Type{Tensor{4, 1}}){T <: Real}(t::SVector{1, T}) = Tensor{4, 1, T, 1}(t)
(::Type{Tensor{4, 2}}){T <: Real}(t::SVector{16, T}) = Tensor{4, 2, T, 16}(t)
(::Type{Tensor{4, 3}}){T <: Real}(t::SVector{81, T}) = Tensor{4, 3, T, 81}(t)

###############
# Typealiases #
###############

typealias Vec{dim, T, M} Tensor{1, dim, T, dim}

typealias AllTensors{dim, T} Union{SymmetricTensor{2, dim, T}, Tensor{2, dim, T},
                                   SymmetricTensor{4, dim, T}, Tensor{4, dim, T},
                                   Vec{dim, T}}


typealias SecondOrderTensor{dim, T} Union{SymmetricTensor{2, dim, T}, Tensor{2, dim, T}}
typealias FourthOrderTensor{dim, T} Union{SymmetricTensor{4, dim, T}, Tensor{4, dim, T}}
typealias SymmetricTensors{dim, T} Union{SymmetricTensor{2, dim, T}, SymmetricTensor{4, dim, T}}
typealias Tensors{dim, T} Union{Tensor{2, dim, T}, Tensor{4, dim, T}, Vec{dim, T}}

include("utilities.jl")
include("tuple_utils.jl")
include("symmetric_tuple_linalg.jl")

include("indexing.jl")
include("promotion_conversion.jl")
include("tensor_ops.jl")
include("tensor_ops_errors.jl")
include("symmetric_ops.jl")


##############################
# Utility/Accessor Functions #
##############################

get_data(t::AbstractTensor) = t.data.data
tovector(t::AbstractTensor) = t.data

function tomatrix{dim}(t::Tensor{4, dim})
    N = n_components(Tensor{2,dim})
    return SMatrix{N, N}(tovector(t))
end

function tomatrix{dim}(t::Tensor{2, dim})
    N = n_components(Tensor{1,dim})
    return SMatrix{N, N}(tovector(t))
end

function n_independent_components(dim, issym)
    dim == 1 && return 1
    if issym
        dim == 2 && return 3
        dim == 3 && return 6
    else
        dim == 2 && return 4
        dim == 3 && return 9
    end
    return -1
end

@pure n_components{dim}(::Type{SymmetricTensor{2, dim}}) = dim*dim - div((dim-1)*dim, 2)
@pure function n_components{dim}(::Type{SymmetricTensor{4, dim}})
    n = n_components(SymmetricTensor{2, dim})
    return n*n
end

@pure n_components{order, dim}(::Type{Tensor{order, dim}}) = dim^order

@pure get_main_type{order, dim, T, M}(::Type{SymmetricTensor{order, dim, T, M}}) = SymmetricTensor
@pure get_main_type{order, dim, T, M}(::Type{Tensor{order, dim, T, M}}) = Tensor
@pure get_main_type{order, dim, T}(::Type{SymmetricTensor{order, dim, T}}) = SymmetricTensor
@pure get_main_type{order, dim, T}(::Type{Tensor{order, dim, T}}) = Tensor
@pure get_main_type{order, dim}(::Type{SymmetricTensor{order, dim}}) = SymmetricTensor
@pure get_main_type{order, dim}(::Type{Tensor{order, dim}}) = Tensor

@pure get_base{order, dim, T, M}(::Type{SymmetricTensor{order, dim, T, M}}) = SymmetricTensor{order, dim}
@pure get_base{order, dim, T, M}(::Type{Tensor{order, dim, T, M}}) = Tensor{order, dim}

@pure get_lower_order_tensor{dim, T, M}(S::Type{SymmetricTensor{2, dim, T, M}}) =  SymmetricTensor{2, dim}
@pure get_lower_order_tensor{dim, T, M}(S::Type{Tensor{2, dim, T, M}}) = Tensor{2, dim}
@pure get_lower_order_tensor{dim, T, M}(::Type{SymmetricTensor{4, dim, T, M}}) = SymmetricTensor{2, dim}
@pure get_lower_order_tensor{dim, T, M}(::Type{Tensor{4, dim, T, M}}) = Tensor{2, dim}


############################
# Abstract Array interface #
############################

Base.linearindexing{T <: SymmetricTensor}(::Type{T}) = Base.LinearSlow()
Base.linearindexing{T <: Tensor}(::Type{T}) = Base.LinearFast()

get_type{X}(::Type{Type{X}}) = X


########
# Size #
########

Base.size{dim}(::Vec{dim}) = (dim,)
Base.size{dim}(::SecondOrderTensor{dim}) = (dim, dim)
Base.size{dim}(::FourthOrderTensor{dim}) = (dim, dim, dim, dim)

Base.similar(t::AbstractTensor) = t

# Ambiguity fix
Base.fill(t::AbstractTensor, v::Integer)  = one(typeof(t)) * v
Base.fill(t::AbstractTensor, v::Number) = one(typeof(t)) * v


#########################
# Internal constructors #
#########################

@inline function (Tt::Union{Type{Tensor{order, dim, T, M}}, Type{SymmetricTensor{order, dim, T, M}}}){order, dim, T, M}(data)
    get_base(Tt)(data)
end

@inline function (Tt::Type{Vec{dim}}){dim}(data)
    Tensor{1, dim}(data)
end

# These are some kinda ugly stuff to create different type of constructors.
@generated function (Tt::Type{Tensor{order, dim}}){order, dim}(data::Union{AbstractArray, Tuple})
    # Check for valid orders
    n = n_components(Tensor{order,dim})
    if !(order in (1,2,4))
        throw(ArgumentError("Tensor only supported for order 1, 2, 4"))
    end
    return quote
        if length(data) != $n
            throw(ArgumentError("Wrong number of tuple elements, expected $($n), got $(length(data))"))
        end
        Tensor{order, dim}(to_tuple(NTuple{$n}, data))
    end
end

# These are some kinda ugly stuff to create different type of constructors.
@generated function (Tt::Type{SymmetricTensor{order, dim}}){order, dim}(data::Union{AbstractArray, Tuple})
    n = n_components(Tensor{order,dim})
    m = n_components(SymmetricTensor{order,dim})
    if !(order in (2,4))
        throw(ArgumentError("SymmetricTensor only supported for order 2, 4"))
    end
    return quote
        if length(data) != $n && length(data) != $m
            throw(ArgumentError("Wrong number of tuple elements, expected $($n) or $($m), got $(length(data))"))
        end
        if length(data) == $m
            return SymmetricTensor{order, dim, eltype(data), $m}(to_tuple(NTuple{$m}, data))
        end
        S = Tensor{order, dim}(to_tuple(NTuple{$n}, data))
        return convert(SymmetricTensor{order, dim}, S)
    end
end

@generated function (Tt::Union{Type{Tensor{order, dim}}, Type{SymmetricTensor{order, dim}}}){order, dim}(f::Function)
    # Check for valid orders
    if !(order in (1,2,4))
        throw(ArgumentError("Only tensors of order 1, 2, 4 supported"))
    end

    # Storage format is of rank 1 for vectors and order / 2 for other tensors
    if order == 1 && Tt <: SymmetricTensor
        throw(ArgumentError("SymmetricTensor only supported for order 2, 4"))
    end

    n = n_components(get_type(Tt))

    # Validate that the input array has the correct number of elements.
    if order == 1
        exp = tensor_create(get_main_type(get_type(Tt)){order, dim}, (i) -> :(f($i)))
    elseif order == 2
        exp = tensor_create(get_main_type(get_type(Tt)){order, dim}, (i,j) -> :(f($i, $j)))
    elseif order == 4
        exp = tensor_create(get_main_type(get_type(Tt)){order, dim}, (i,j,k,l) -> :(f($i, $j, $k, $l)))
    end

    return :(get_main_type(Tt){order, dim}($exp))
end

function (Tt::Union{Type{Tensor{order, dim, T}}, Type{SymmetricTensor{order, dim, T}}}){order, dim, T}(f_or_data)
    t1 = get_main_type(Tt){order, dim}(f_or_data)
    return convert(get_main_type(Tt){order, dim}, t1)
end


###############
# Simple Math #
###############

# Num, tensor. *, /
for TensorType in (SymmetricTensor, Tensor)
    @eval begin
        @inline Base.:*{order, dim, T, N}(n::Number, t::$TensorType{order, dim, T, N}) = $TensorType{order, dim}(n * tovector(t))
        @inline Base.:*{order, dim, T, N}(t::$TensorType{order, dim, T, N}, n::Number) = $TensorType{order, dim}(tovector(t) * n)
        @inline Base.:/{order, dim, T, N}(t::$TensorType{order, dim, T, N}, n::Number) = $TensorType{order, dim}(tovector(t) / n)

        # Unary -, +
        @inline Base.:-{order, dim, T, N}(t::$TensorType{order, dim, T, N}) = $TensorType{order, dim}(-tovector(t))
        @inline Base.:+{order, dim, T, N}(t::$TensorType{order, dim, T, N}) = $TensorType{order, dim}(+tovector(t))
    end
end

# Binary, +, -, -*, ./
for (op) in (:-, :+, :.*, :./, :.-, :.+)
    for TensorType in (SymmetricTensor, Tensor)
        @eval begin
            @inline function Base.$op{order, dim, T1, T2, N}(t1::$TensorType{order, dim, T1, N}, t2::$TensorType{order, dim, T2, N})
                $TensorType{order, dim}($op(tovector(t1), tovector(t2)))
            end
        end
    end
    @eval begin
        Base.$op{order, dim}(t1::AbstractTensor{order, dim}, t2::AbstractTensor{order, dim}) = Base.$op(promote(t1, t2)...)
    end
end

# zero, rand
for op in (:zero, :rand)
    for TensorType in (SymmetricTensor, Tensor)
        @eval begin
            @inline function Base.$op{order, dim, T}(Tt::Type{$TensorType{order, dim, T}})
                N = n_components($TensorType{order, dim})
                return $TensorType{order, dim}($op(SVector{N, T}))
            end
        end
    end
end


for op in (:zero, :rand, :one)
    for TensorType in (SymmetricTensor, Tensor)
        @eval begin
            @inline Base.$op{order, dim}(Tt::Type{$TensorType{order, dim}}) = Base.$op($TensorType{order, dim, Float64})

            @inline function Base.$op{order, dim, T, M}(Tt::Type{$(TensorType){order, dim, T, M}})
                $op($(TensorType){order, dim})
            end
        end
    end

    # Special for Vecs
    @eval begin
        @inline function Base.$op{dim}(Tt::Type{Vec{dim}})
            $op(Vec{dim, Float64})
        end

        @inline function Base.$op(t::AllTensors)
            $op(typeof(t))
        end
    end
end

##########################
# Zero, one, rand, diagm #
##########################

for TensorType in (SymmetricTensor, Tensor)
    @eval begin
        @generated function Base.diagm{order, dim, T}(Tt::Type{$(TensorType){order, dim}}, v::AbstractVector{T})
            N = n_components($(TensorType){order, dim})
            if order == 1
                f = (i) -> :(v[$i])
            elseif order == 2
                f = (i,j) -> i == j ? :(v[$i]) : :($(zero(T)))
            elseif order == 4
                f = (i,j,k,l) -> i == k && j == l ? :(v[$i]) : :($(zero(T)))
            end
            exp = tensor_create(get_type(Tt),f)
            return quote
                $(Expr(:meta, :inline))
                @inbounds t = $exp
                $($TensorType){order, dim}(t)
            end
        end

        @generated function Base.diagm{order, dim, T}(Tt::Type{$(TensorType){order, dim}}, v::T)
            N = n_components($(TensorType){order, dim})
            if order == 1
                f = (i) -> :(v)
            elseif order == 2
                f = (i,j) -> i == j ? :(v) : :($(zero(T)))
            elseif order == 4
                f = (i,j,k,l) -> i == k && j == l ? :(v) : :($(zero(T)))
            end
            exp = tensor_create(get_type(Tt),f)
            return quote
                $(Expr(:meta, :inline))
                $($TensorType){order, dim}($exp)
            end
        end

        function Base.one{order, dim, T}(Tt::Type{$(TensorType){order, dim, T}})
            Base.diagm($(TensorType){order, dim}, one(T))
        end

    end
end

end # module
