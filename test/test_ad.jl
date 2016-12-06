
for dim in 1:3
    for T in (Float32, Float64)
        A = rand(Tensor{2, dim, T})
        A_sym = rand(SymmetricTensor{2, dim, T})
        v = rand(Vec{dim, T})

        # Gradient of scalars
        @test ContMechTensors.gradient(norm, v) ≈ v / norm(v)
        @test ContMechTensors.gradient(v -> 3*v, v) ≈ diagm(Tensor{2, dim}, 3.0)
        # https://en.wikipedia.org/wiki/Tensor_derivative_(continuum_mechanics)#Derivatives_of_the_invariants_of_a_second-order_tensor
        I1, DI1 = A -> trace(A), A -> one(A)
        I2, DI2 = A -> 1/2 * (trace(A)^2 - trace(A⋅A)), A -> I1(A) * one(A) - A'
        I3, DI3 = A -> det(A), A -> det(A) * inv(A)'

        @test ContMechTensors.gradient(I1, A) ≈ DI1(A)
        @test ContMechTensors.gradient(I2, A) ≈ DI2(A)
        @test ContMechTensors.gradient(I3, A) ≈ DI3(A)
        @test ContMechTensors.gradient(I1, A_sym) ≈ DI1(A_sym)
        @test ContMechTensors.gradient(I2, A_sym) ≈ DI2(A_sym)
        @test ContMechTensors.gradient(I3, A_sym) ≈ DI3(A_sym)
    end
end