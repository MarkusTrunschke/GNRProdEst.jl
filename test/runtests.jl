using GNRProdEst
using Test
using CSV
using DataFrames
using Statistics

@testset "GNRProdEst.jl" begin

        # Get data and prepare for tests
        # data = CSV.read("C:/Users/marku/Documents/GNRProdEst/Other Programs/R-version/cd_data_500.csv", DataFrame)
        data = CSV.read("C:/Users/marku/Documents/GNRProdEst/Other Programs/GNR/Cleaned version of Table_1/cd_data_500.csv", DataFrame)
        data[!, :share_m_y] = select(data, :i_level)[:,1] ./ select(data, :yg_level)[:,1] # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...
        data[!, :ln_share_m_y] = log.(select(data, :i_level)[:,1] ./ select(data, :yg_level)[:,1]) # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...
        data.constant = ones(size(data)[1])

        data_mis = allowmissing(data)
        data_mis[1,5] = missing
        opts = Dict( "fs_series_order" => 2,
                     "fs_method" => "NLLS", # OLS is experimental. It gives (almost) the same results but the residual is different
                     "fs_print_starting_values" => true,
                     "fs_print_results" => true)


        fes_res, ses_res = GNRProdEst.GNRProd(data = data_mis, output = :yg, flex_input = :i, fixed_inputs = :k, id = :id, time=:time, opts = opts)

        γ_dash_GNR = [0.6518793, 0.0034517, -0.0011976, 0.0013776, -0.0046622, -0.0000291]
        γ_dash_R = [6.519479e-01, 3.464846e-03, -1.196053e-03, 1.398264e-03, -4.688621e-03, -5.307226e-05]

        result_df = select(result_df, Not(:flex_elas2))
        γ, γ_flex, E = GNRProdEst.fes_predictions!(data = result_df, ln_share_flex_y_var = :ln_share_flex_y, flex_input = :i, input_var_symbols=fes_res["taylor_series"], γ_dash = γ_dash_R, output = :yg)



        data = CSV.read("C:/Users/marku/Documents/GNRProdEst/Other Programs/R-version/Columbia_311.csv", DataFrame)
        data.constant = ones(size(data)[1])

        # GNRProdEst.GNRFirstStage(est_df = data, output = :Y, flex_inputs = :M, fixed_inputs = :K, share = :ln_share_m_y, all_input_symbols = [:K, :M], opts = opts)
        flex_elas, ln_int_G_I, share, ϵ = GNRProdEst.GNRFirstStage(est_df = data, output = :yg, flex_input = :i, fixed_inputs = :k, share = :si, all_input_symbols = [:k, :i], id = :id, time = :time, opts = opts);
        # flex_elas, ln_int_G_I, share, ϵ = GNRProdEst.GNRFirstStage(est_df = data, output = :RGO, flex_input = :RI, fixed_inputs = [:L, :K], share = :share, all_input_symbols = [:L, :K, :RI], opts = opts);

        # Test data preperation (drop missings and add column)
        @test GNRProdEst.prep_data(data_mis; output = :Y, flex_inputs = :M, fixed_inputs = :K, intermediate_input = :M)[1] == dropmissing(data_mis[:, [:Y, :M, :K, :constant, :ln_share_m_y]])
        @test GNRProdEst.GNRFirstStage(est_df = data, output = :Y, flex_inputs = :M, fixed_inputs = :K, share = :ln_share_m_y, all_input_symbols = [:K, :M], opts = opts) ≈ [ -0.0010269005340738705, 5.573265041379824e-7, -8.676082093980585e-8, -0.0001222493225217935, 3.6339341714268326e-8]

        GNRProd(;data::DataFrame, output::Symbol, flex_inputs::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, opts::Dict)
    
end




opts = Dict( "fes_series_order" => 2,
             "fes_method" => "NLLS", # OLS is experimental. It gives (almost) the same results but the residual is different
             "fes_print_starting_values" => true,
             "fes_print_results" => true,
             "ses_series_order" => 2)

GNRProdEst.GNRSecondStage(est_df = data, fixed_inputs = [:k], called_from_GNRProd = true, fes_returns = fes_res, opts = opts)




new_df = GNRProdEst.panel_lag(;data = data, id = :id, time = :time, variable = :k, lag_prefix = "lag_", lags = 1, drop_missings = false)
