using GNRProdEst
using Test
using CSV
using DataFrames

@testset "GNRProdEst.jl" begin

        # Get data and prepare for tests
        data = CSV.read("C:/Users/marku/Documents/GNRProdEst/Other Programs/R-version/cd_data_500.csv", DataFrame)
        data_mis = allowmissing(data)
        data_mis[1,5] = missing
        data_mis[!, :ln_share_m_y] = log.(select(data_mis, :M)[:,1] ./ select(data_mis, :Y)[:,1]) # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...

        # Test data preperation (drop missings and add column)
        @test GNRProdEst.prep_data(data_mis; output = :Y, flex_inputs = :M, fixed_inputs = :K, intermediate_input = :M) == dropmissing!(data_mis[:, [:Y, :M, :K, :ln_share_m_y]])


end
