using GNRProdEst
using Test
using CSV
using DataFrames

@testset "GNRProdEst.jl" begin
        data = CSV.read("C:/Users/marku/Documents/GNRProdEst/Other Programs/R-version/cd_data.csv", DataFrame)
        @test sin(-θ) ≈ -sin(θ)
        @test cos(-θ) ≈ cos(θ)
        @test sin(2θ) ≈ 2*sin(θ)*cos(θ)
        @test cos(2θ) ≈ cos(θ)^2 - sin(θ)^2
end
