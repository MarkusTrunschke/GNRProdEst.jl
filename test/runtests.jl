using GNRProdEst
using Test
using CSV
using DataFrames

@testset "GNRProdEst.jl" begin

        # Get data and prepare for tests
        # data = CSV.read("C:/Users/marku/Documents/GNRProdEst/Other Programs/R-version/cd_data_500.csv", DataFrame)
        data = CSV.read("C:/Users/marku/Documents/GNRProdEst/Other Programs/GNR/Cleaned version of Table_1/cd_data_500.csv", DataFrame)
        data[!, :ln_share_m_y] = log.(select(data_mis, :M)[:,1] ./ select(data_mis, :Y)[:,1]) # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...
        data.constant = ones(size(data)[1])

        data_mis = allowmissing(data)
        data_mis[1,5] = missing
        opts = Dict( "fs_series_order" => 2,
                     "fs_method" => "OLS",
                     "fs_print_starting_values" => false,
                     "fs_print_results" => true)

        # GNRProdEst.GNRFirstStage(est_df = data, output = :Y, flex_inputs = :M, fixed_inputs = :K, share = :ln_share_m_y, all_input_symbols = [:K, :M], opts = opts)
        OLS = GNRProdEst.GNRFirstStage(est_df = data, output = :yg, flex_inputs = :i, fixed_inputs = :k, share = :si, all_input_symbols = [:k, :i], opts = opts)
        
        # Test data preperation (drop missings and add column)
        @test GNRProdEst.prep_data(data_mis; output = :Y, flex_inputs = :M, fixed_inputs = :K, intermediate_input = :M)[1] == dropmissing(data_mis[:, [:Y, :M, :K, :constant, :ln_share_m_y]])
        @test GNRProdEst.GNRFirstStage(est_df = data, output = :Y, flex_inputs = :M, fixed_inputs = :K, share = :ln_share_m_y, all_input_symbols = [:K, :M], opts = opts) ≈ [ -0.0010269005340738705, 5.573265041379824e-7, -8.676082093980585e-8, -0.0001222493225217935, 3.6339341714268326e-8]

        GNRProd(;data::DataFrame, output::Symbol, flex_inputs::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, opts::Dict)
    
end

