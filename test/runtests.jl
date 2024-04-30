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
                     "fs_method" => "NLLS", # OLS is experimental. It gives (almost) the same results but the residual is different
                     "fs_print_starting_values" => true,
                     "fs_print_results" => true)

        data = CSV.read("C:/Users/marku/Documents/GNRProdEst/Other Programs/R-version/Columbia_311.csv", DataFrame)
        data.constant = ones(size(data)[1])

        gnr_FS <- gnrflex(output = "RGO",
                   fixed = c("L", "K"),
                   flex = "RI",
                   share = "share",
                   id = "id",
                   time = "year",
                   data = Columbia_311,
                   control = list(degree = 2, maxit = 2000))
        # GNRProdEst.GNRFirstStage(est_df = data, output = :Y, flex_inputs = :M, fixed_inputs = :K, share = :ln_share_m_y, all_input_symbols = [:K, :M], opts = opts)
        flex_elas, ln_int_G_I, share, ϵ = GNRProdEst.GNRFirstStage(est_df = data, output = :yg, flex_input = :i, fixed_inputs = :k, share = :si, all_input_symbols = [:k, :i], opts = opts);
        # flex_elas, ln_int_G_I, share, ϵ = GNRProdEst.GNRFirstStage(est_df = data, output = :RGO, flex_input = :RI, fixed_inputs = [:L, :K], share = :share, all_input_symbols = [:L, :K, :RI], opts = opts);

        # Test data preperation (drop missings and add column)
        @test GNRProdEst.prep_data(data_mis; output = :Y, flex_inputs = :M, fixed_inputs = :K, intermediate_input = :M)[1] == dropmissing(data_mis[:, [:Y, :M, :K, :constant, :ln_share_m_y]])
        @test GNRProdEst.GNRFirstStage(est_df = data, output = :Y, flex_inputs = :M, fixed_inputs = :K, share = :ln_share_m_y, all_input_symbols = [:K, :M], opts = opts) ≈ [ -0.0010269005340738705, 5.573265041379824e-7, -8.676082093980585e-8, -0.0001222493225217935, 3.6339341714268326e-8]

        GNRProd(;data::DataFrame, output::Symbol, flex_inputs::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, opts::Dict)
    
end


data = CSV.read("C:/Users/marku/Documents/GNRProdEst/Other Programs/GNR/Cleaned version of Table_1/cd_data_500.csv", DataFrame)

est_df, all_var_symbols, all_input_symbols = GNRProdEst.prep_data(data, output = :yg , flex_input = :i, fixed_inputs = :k)

data[!, :share_flex_y] = data.i_level ./ data.yg_level # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...
data[!, :ln_share_flex_y] = data.i .- data.yg # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...
data.constant = ones(size(data)[1])

coefs = [ 6.519479e-01, -4.688618e-03, -5.307394e-05, 1.398266e-03, 3.464840e-03,  -1.196054e-03]

input_var_symbols = [:k, :k_k, :k_i, :i, :i_i]


flex_elas[1:6]
