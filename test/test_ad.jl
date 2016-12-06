import ContMechTensors.gradient
for dim in 1:3
    for T in (Float32, Float64)
        c = rand(Tensor{2, dim, T})
        v = rand(Vec{dim, T})
        @test gradient(norm, v) ≈ v / norm(v)
        @test gradient(v -> 3*v, v) ≈ diagm(Tensor{2, dim}, 3.0)
    end
end