var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ContMechTensors-1",
    "page": "Home",
    "title": "ContMechTensors",
    "category": "section",
    "text": "Efficient computations with symmetric and unsymmetric tensors in Julia."
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "This Julia package provides fast operations with symmetric/unsymmetric tensors of order 1, 2 and 4. The tensors are stack allocated which means that there is no need to preallocate results of operations and nice infix notation can be used without a performance penalty. For the symmetric tensors, when possible, the symmetry is exploited for better performance."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "ContMechTensors is a registered package and so can be installed viaPkg.add(\"ContMechTensors\")The package has no dependencies other than Julia (0.5 and up) itself."
},

{
    "location": "index.html#Manual-Outline-1",
    "page": "Home",
    "title": "Manual Outline",
    "category": "section",
    "text": "Pages = [\n    \"man/constructing_tensors.md\",\n    \"man/indexing.md\",\n    \"man/binary_operators.md\",\n    \"man/other_operators.md\",\n    \"man/storing_tensors.md\",\n]\nDepth = 1"
},

{
    "location": "index.html#Demos-1",
    "page": "Home",
    "title": "Demos",
    "category": "section",
    "text": "Pages = [\n    \"demos.md\"]\nDepth = 1"
},

{
    "location": "man/constructing_tensors.html#",
    "page": "Constructing tensors",
    "title": "Constructing tensors",
    "category": "page",
    "text": "DocTestSetup = quote\n    srand(1234)\n    using ContMechTensors\nend"
},

{
    "location": "man/constructing_tensors.html#Constructing-tensors-1",
    "page": "Constructing tensors",
    "title": "Constructing tensors",
    "category": "section",
    "text": "Tensors can be created in multiple ways but they usually include Tensor{order, dim} or SymmetricTensor{order, dim}julia> zero(Tensor{1, 2})\n2-element ContMechTensors.Tensor{1,2,Float64,2}:\n 0.0\n 0.0\n\njulia> rand(Tensor{2, 3})\n3×3 ContMechTensors.Tensor{2,3,Float64,9}:\n 0.590845  0.460085  0.200586\n 0.766797  0.794026  0.298614\n 0.566237  0.854147  0.246837\n\njulia> zero(SymmetricTensor{4, 2})\n2×2×2×2 ContMechTensors.SymmetricTensor{4,2,Float64,9}:\n[:, :, 1, 1] =\n 0.0  0.0\n 0.0  0.0\n\n[:, :, 2, 1] =\n 0.0  0.0\n 0.0  0.0\n\n[:, :, 1, 2] =\n 0.0  0.0\n 0.0  0.0\n\n[:, :, 2, 2] =\n 0.0  0.0\n 0.0  0.0\n\njulia> one(SymmetricTensor{2, 2})\n2×2 ContMechTensors.SymmetricTensor{2,2,Float64,3}:\n 1.0  0.0\n 0.0  1.0Tensors can also be created by giving a tuple or an array with the same number of elements as the number of independent indices in the tensor:julia> Tensor{1,2}([1.0,2.0])\n2-element ContMechTensors.Tensor{1,2,Float64,2}:\n 1.0\n 2.0\n\njulia> SymmetricTensor{2,2}((1.0,2.0,3.0))\n2×2 ContMechTensors.SymmetricTensor{2,2,Float64,3}:\n 1.0  2.0\n 2.0  3.0It is also possible to create a tensor by giving a function f(index...) -> v:julia> SymmetricTensor{2,2}((i,j) -> i + j)\n2×2 ContMechTensors.SymmetricTensor{2,2,Int64,3}:\n 2  3\n 3  4A diagonal tensor can be created by either giving a number of a vector on the diagonal:julia> diagm(Tensor{2,2}, 2.0)\n2×2 ContMechTensors.Tensor{2,2,Float64,4}:\n 2.0  0.0\n 0.0  2.0\n\njulia> diagm(SymmetricTensor{2,3}, [1.0, 2.0, 3.0])\n3×3 ContMechTensors.SymmetricTensor{2,3,Float64,6}:\n 1.0  0.0  0.0\n 0.0  2.0  0.0\n 0.0  0.0  3.0"
},

{
    "location": "man/indexing.html#",
    "page": "Indexing",
    "title": "Indexing",
    "category": "page",
    "text": "DocTestSetup = quote\n    srand(1234)\n    using ContMechTensors\nend"
},

{
    "location": "man/indexing.html#Indexing-1",
    "page": "Indexing",
    "title": "Indexing",
    "category": "section",
    "text": "Indexing into a (Symmetric)Tensor{dim, order} is performed like for an Array of dimension order.julia> A = rand(Tensor{2, 2});\n\njulia> A[1, 2]\n0.5662374165061859\n\njulia> B = rand(SymmetricTensor{4, 2});\n\njulia> B[1, 2, 1, 2]\n0.24683718661000897In order to set an index the function setindex(t, value, index...) is used. This returns a new tensor with the modified index. Explicitly setting indices is not recommended in performance critical code since it will invoke dynamic dispatch. It is provided as a means of convenience when working in for example the REPL.julia> a = rand(Vec{2});\n\njulia> setindex(a, 1.337, 2)\n2-element ContMechTensors.Tensor{1,2,Float64,2}:\n 0.590845\n 1.337"
},

{
    "location": "man/binary_operators.html#",
    "page": "Binary Operations",
    "title": "Binary Operations",
    "category": "page",
    "text": "DocTestSetup = quote\n    srand(1234)\n    using ContMechTensors\nend"
},

{
    "location": "man/binary_operators.html#Binary-Operations-1",
    "page": "Binary Operations",
    "title": "Binary Operations",
    "category": "section",
    "text": ""
},

{
    "location": "man/binary_operators.html#Single-contraction-(dot-product)-1",
    "page": "Binary Operations",
    "title": "Single contraction (dot product)",
    "category": "section",
    "text": "Single contractions or scalar products of a tensor with order n and a tensor with order m gives a tensor with order m + n - 2. The symbol ⋅, written \\cdot, is overloaded for single contraction.julia> A = rand(Tensor{2, 2})\n2×2 ContMechTensors.Tensor{2,2,Float64,4}:\n 0.590845  0.566237\n 0.766797  0.460085\n\njulia> B = rand(Tensor{1, 2})\n2-element ContMechTensors.Tensor{1,2,Float64,2}:\n 0.794026\n 0.854147\n\njulia> dot(A, B)\n2-element ContMechTensors.Tensor{1,2,Float64,2}:\n 0.952796\n 1.00184\n\njulia> A ⋅ B\n2-element ContMechTensors.Tensor{1,2,Float64,2}:\n 0.952796\n 1.00184"
},

{
    "location": "man/binary_operators.html#Double-contraction-1",
    "page": "Binary Operations",
    "title": "Double contraction",
    "category": "section",
    "text": "Double contractions contracts the two most inner \"legs\" of the tensors. The result of a double contraction between a tensor of order n and a tensor with order m gives a tensor with order m + n - 4. The symbol ⊡, written \\boxdot, is overloaded for double contraction. The reason : is not used is because it does not have the same precedence as multiplication.julia> A = rand(SymmetricTensor{2, 2});\n\njulia> B = rand(SymmetricTensor{2, 2});\n\njulia> dcontract(A,B)\n1.9732018397544984\n\njulia> A ⊡ B\n1.9732018397544984"
},

{
    "location": "man/binary_operators.html#Tensor-product-(open-product)-1",
    "page": "Binary Operations",
    "title": "Tensor product (open product)",
    "category": "section",
    "text": "Tensor products or open product of a tensor with order n and a tensor with order m gives a tensor with order m + n. The symbol ⊗, written \\otimes, is overloaded for tensor products.julia> A = rand(SymmetricTensor{2, 2});\n\njulia> B = rand(SymmetricTensor{2, 2});\n\njulia> A ⊗ B\n2×2×2×2 ContMechTensors.SymmetricTensor{4,2,Float64,9}:\n[:, :, 1, 1] =\n 0.271839  0.352792\n 0.352792  0.260518\n\n[:, :, 2, 1] =\n 0.469146  0.608857\n 0.608857  0.449607\n\n[:, :, 1, 2] =\n 0.469146  0.608857\n 0.608857  0.449607\n\n[:, :, 2, 2] =\n 0.504668  0.654957\n 0.654957  0.48365"
},

{
    "location": "man/other_operators.html#",
    "page": "Other operators",
    "title": "Other operators",
    "category": "page",
    "text": ""
},

{
    "location": "man/other_operators.html#Other-operators-1",
    "page": "Other operators",
    "title": "Other operators",
    "category": "section",
    "text": "For vectors (first order tensors): normFor second order tensors: norm, trace (vol), det, inv, transpose, symmetric, skew, eig, mean defined as trace(s) / 3, and dev defined as s - mean(s) * I.For fourth order tensors: norm, trace, symmetric (same as minorsymmetric), majorsymmetric, transpose (same as minortranspose), majortranspose, permute_indexThere is also a few special functions that can be convenient:For computing F' ⋅ F between two general second order tensors there is tdot(F) which returns a SymmetricTensor.\nFor computing a_k cdot mathbfC_ikjl cdot b_l for two vectors a and b and a fourth order symmetric tensor mathbfC there is dotdot(a, C, b). This function is useful because it is the expression for the tangent matrix in continuum mechanics when the displacements are approximated by scalar base functions."
},

{
    "location": "man/storing_tensors.html#",
    "page": "Storing tensors",
    "title": "Storing tensors",
    "category": "page",
    "text": ""
},

{
    "location": "man/storing_tensors.html#Storing-tensors-1",
    "page": "Storing tensors",
    "title": "Storing tensors",
    "category": "section",
    "text": "Even though a user mostly deals with the Tensor{order, dim, T} parameters, the full parameter list for a tensor is actually Tensor{order, dim, T, N} where N is the number of independent elements in the tensor. The reason for this is that the internal storage is a NTuple{N, T}. In order to get good performance when storing tensors in other types it is importatant that the container type is also parametrized on N. For example, when storing one symmetric second order tensor and one unsymmetric tensor, this is the preferred way:immutable Container{dim, T, N, M}\n    sym_tens::SymmetricTensor{2, dim, T, N}\n    tens::Tensor{2, dim, T, M}\nendLeaving out the M and N would lead to bad performance."
},

{
    "location": "demos.html#",
    "page": "Demos",
    "title": "Demos",
    "category": "page",
    "text": ""
},

{
    "location": "demos.html#Demos-1",
    "page": "Demos",
    "title": "Demos",
    "category": "section",
    "text": "This section contain a few demos of applying ContMechTensors to continuum mechanics."
},

{
    "location": "demos.html#Creating-the-linear-elasticity-tensor-1",
    "page": "Demos",
    "title": "Creating the linear elasticity tensor",
    "category": "section",
    "text": "The linear elasticity tensor mathbfC can be defined from the Lame parameters lambda and mu by the expression$ \\mathbf{C}_{ijkl} = \\lambda \\delta_{ij}\\delta_{kl} + \\mu(\\delta_{ij}\\delta_{jl} + \\delta_{il}\\delta_{jk}),$where delta_ij = 1 if i = j otherwise 0. It can also be computed in terms of the Young's modulus E and Poisson's ratio nu by the conversion formulas lambda = Enu  (1 + nu)(1 - 2nu) and mu = E  2(1 + nu)The code below creates the elasticity tensor for given parameters E and nu and dimension textttdim. Note the similarity between the mathematical formula and the code.using ContMechTensors\nE = 200e9\nν = 0.3\ndim = 2\nλ = E*ν / ((1 + ν) * (1 - 2ν))\nμ = E / (2(1 + ν))\nδ(i,j) = i == j ? 1.0 : 0.0\nf = (i,j,k,l) -> λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))\n\nC = SymmetricTensor{4, dim}(f)"
},

{
    "location": "demos.html#Nonlinear-elasticity-material-1",
    "page": "Demos",
    "title": "Nonlinear elasticity material",
    "category": "section",
    "text": "For a deformation gradient mathbfF = mathbfI + nabla otimes mathbfu, where mathbfu is the deformation from the reference to the current configuration, the right Cauchy-Green deformation tensor is defined by mathbfC = mathbfF^T cdot mathbfF. The Second Piola Krichoff stress tensor mathbfS is derived from the Helmholtz free energy Psi by the relation mathbfS = 2 fracpartial Psipartial mathbfC.We can define the energy for a material with thePsi(mathbfC) = 12 mu (mathrmtr(hatmathbfC) - 3) + K_b(J-1)^2where hatmathbfC = mathrmdet(mathbfC)^-13 mathbfC and J = det(mathbfF) = sqrtdet(mathbfC) and the shear and bulk modulus are given by mu and K_b respectively.This free energy function can be implemented as:function Ψ(C, μ, Kb)\n    detC = det(C)\n    J = sqrt(detC)\n    Ĉ = detC^(-1/3)*C\n    return 1/2*(μ * (trace(Ĉ)- 3) + Kb*(J-1)^2)\nendThe analytical expression for the Second Piola Kirchoff tensor is$ \\mathbf{S} = \\mu \\det(\\mathbf{C})^{-1/3}(\\mathbf{I} - 1/3 \\mathrm{tr}(\\mathbf{C})\\mathbf{C}^{-1}) + K_b(J-1)J\\mathbf{C}^{-1} $which can be implemented by the functionfunction S(C, μ, Kb)\n    I = one(C)\n    J = sqrt(det(C))\n    invC = inv(C)\n    return μ * det(C)^(-1/3)*(I - 1/3*trace(C)*invC) + Kb*(J-1)*J*invC\nend"
},

{
    "location": "demos.html#Automatic-differentiation-1",
    "page": "Demos",
    "title": "Automatic differentiation",
    "category": "section",
    "text": "For some material models it can be cumbersome to compute the analytical expression for the Second Piola Kirchoff tensor. We can then use Automatic Differentiation (AD) to compute it. Here, the AD package (ForwardDiff.jl)[https://github.com/JuliaDiff/ForwardDiff.jl] is used. Unfortunately we have to here do a bit of juggling between tensors and standard Julia Arrays due to ForwardDiff expecting the input to be of Array type.using ForwardDiff\n\nfunction S_AD{dim}(C::SymmetricTensor{2,dim}, μ, Kb)\n    Ψvec = Cvec -> Ψ(SymmetricTensor{2,dim}(Cvec), μ, Kb)\n    ∂Ψ∂C = C -> symmetric(Tensor{2,dim}(ForwardDiff.gradient(Ψvec, vec(C))))\n    return 2 * ∂Ψ∂C(C)\nendWe can compare the results from the analytical and AD functions and they are obviously equal:DocTestSetup = quote\n    srand(1234)\n    using ContMechTensors\n    E = 200e9\n    ν = 0.3\n    dim = 2\n    λ = E*ν / ((1 + ν) * (1 - 2ν))\n    μ = E / (2(1 + ν))\n    δ(i,j) = i == j ? 1.0 : 0.0\n    f = (i,j,k,l) -> λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))\n\n    C = SymmetricTensor{4, dim}(f)\n\n    function Ψ(C, μ, Kb)\n        detC = det(C)\n        J = sqrt(detC)\n        Ĉ = detC^(-1/3)*C\n        return 1/2*(μ * (trace(Ĉ)- 3) + Kb*(J-1)^2)\n    end\n\n    function S(C, μ, Kb)\n        I = one(C)\n        J = sqrt(det(C))\n        invC = inv(C)\n        return μ * det(C)^(-1/3)*(I - 1/3*trace(C)*invC) + Kb*(J-1)*J*invC\n    end\n\n    using ForwardDiff\n\n    function S_AD{dim}(C::SymmetricTensor{2,dim}, μ, Kb)\n        Ψvec = Cvec -> Ψ(SymmetricTensor{2,dim}(Cvec), μ, Kb)\n        ∂Ψ∂C = C -> symmetric(Tensor{2,dim}(ForwardDiff.gradient(Ψvec, vec(C))))\n        return 2 * ∂Ψ∂C(C)\n    end\n\nendjulia> μ = 1e10;\n\njulia> Kb = 1.66e11;\n\njulia> F = one(Tensor{2,3}) + rand(Tensor{2,3});\n\njulia> C = tdot(F);\n\njulia> S_AD(C, μ, Kb)\n3×3 ContMechTensors.SymmetricTensor{2,3,Float64,6}:\n  4.30534e11  -2.30282e11  -8.52861e10\n -2.30282e11   4.38793e11  -2.64481e11\n -8.52861e10  -2.64481e11   7.85515e11\n\njulia> S(C, μ, Kb)\n3×3 ContMechTensors.SymmetricTensor{2,3,Float64,6}:\n  4.30534e11  -2.30282e11  -8.52861e10\n -2.30282e11   4.38793e11  -2.64481e11\n -8.52861e10  -2.64481e11   7.85515e11"
},

]}
