import ContMechTensors.gradient
for dim in 1:3
    for T in (Float32, Float64)
        v = rand(Vec{dim, T})
        @test gradient(norm, v) ≈ v / norm(v)
    end
end