using GNRProdEst
using Test
using CSV
using DataFrames
using Statistics
using Profile, ProfileView
using Optim

@testset "GNRProdEst.jl" begin

end




opts = Dict( "fes_method" => "NLLS", # OLS is experimental. It gives (almost) the same results but the residual is slightly different. However, it does not seem to matter too much when you check the second stage results
             "fes_print_starting_values" => true,
             "fes_print_results" => true,
             "ses_optimizer_options" =>  Optim.Options(
                                        f_tol = 1e-9,
                                        x_tol = 1e-2,
                                        g_tol = 1e-10,),
             "ses_print_starting_values" => true,
             "ses_print_results" => true)


# # Get fs results from R for comparison
data_R = CSV.read("/Users/markus_trunschke/Documents/GNRProdEst/Other Programs/GNR/Cleaned version of Table_1/cd_data_500_w_fs.csv", DataFrame)
data_R.weird_Y = data_R.big_Y

data_R.constant = ones(length(data_R.id))
data_R.k2 = data_R.k .* rand(Float64, length(data_R.k))

fes_res, est_df = GNRProdEst.gnrfirststage(est_df = data_R, output = :yg, flexible_input = :i, fixed_inputs = [:k], ln_share_flex_y_var = :si, share_degree = 2, opts = opts);
ses_res, est_df = GNRProdEst.gnrsecondstage(est_df = est_df, id = :id, time = :time, fixed_inputs = :k, flexible_input = :i, fes_returns = fes_res, int_const_series_degree = 2, lm_tfp_degree = 2, opts = opts);

fes_res = GNRProdEst.gnrfirststage!(est_df = data_R, output = :yg, flexible_input = :i, fixed_inputs = [:k], ln_share_flex_y_var = :si, all_input_symbols = [:k, :i], share_degree = 2, opts = opts);
ses_res = GNRProdEst.gnrsecondstage!(est_df = data_R, id = :id, time = :time, fixed_inputs = :k, flexible_input = :i, fes_returns = fes_res, int_const_series_degree = 2, lm_tfp_degree = 2, opts = opts);


gnr_res = GNRProdEst.gnrprodest!(data = data_R, 
                                output = :yg, 
                                flexible_input = :i, 
                                fixed_inputs = :k, 
                                ln_share_flex_y_var = :si, 
                                id = :id, 
                                time = :time,
                                share_degree = 3,
                                int_const_series_degree = 3, 
                                lm_tfp_degree = 3,
                                opts = opts);










# Test with R-version's Colombian test data (likely industry 311)
colombianR_data = CSV.read("/Users/markus_trunschke/Documents/GNRProdEst/Other Programs/R-version/colombian.csv", DataFrame)
colombianR_data.constant = ones(length(colombianR_data.RGO))


opts = Dict("fes_series_degree" => 3,
            "fes_method" => "OLS", # OLS is experimental. It gives (almost) the same results but the residual is slightly different. However, it does not seem to matter too much when you check the second stage results
            "fes_print_starting_values" => true,
            "fes_print_results" => true)


opts = Dict( "fes_method" => "NLLS", # OLS is experimental. It gives (almost) the same results but the residual is slightly different. However, it does not seem to matter too much when you check the second stage results
            "fes_print_starting_values" => true,
            "fes_print_results" => true,
            "fes_optimizer_options" =>  Optim.Options(
                                       f_tol = 1e-9,
                                       x_tol = 1e-8,
                                       g_tol = 1e-10,),
            "ses_optimizer" => IPNewton(),
            "ses_optimizer_options" =>  Optim.Options(
                                       f_tol = 1e-9,
                                       x_tol = 1e-2,
                                       g_tol = 1e-10,),
            "ses_print_starting_values" => true,
            "ses_print_results" => true)



opts = Dict("fes_series_degree" => 3,
            "lm_tfp_degree" => 3,
            "int_const_series_degree" => 3,
            "fes_method" => "NLLS", # OLS is experimental. It gives (almost) the same results but the residual is slightly different. However, it does not seem to matter too much when you check the second stage results
            "fes_print_starting_values" => true,
            "fes_print_results" => true,
            "ses_print_starting_values" => true,
            "ses_print_results" => true)

fes_res, est_df = GNRProdEst.gnrfirststage(est_df = colombianR_data, output = :RGO, flexible_input = :RI, fixed_inputs = [:L :K], ln_share_flex_y_var = :share, share_degree = 3, opts = opts);
ses_res = GNRProdEst.gnrsecondstage(est_df = est_df, id = :id, time = :year, fixed_inputs = [:L :K], flexible_input = :RI, fes_returns = fes_res, int_const_series_degree = 3, lm_tfp_degree = 3, opts = opts);


fes_res = GNRProdEst.gnrfirststage!(est_df = colombianR_data, output = :RGO, flexible_input = :RI, fixed_inputs = [:L :K], ln_share_flex_y_var = :share, share_degree = 3, opts = opts);
ses_res = GNRProdEst.gnrsecondstage!(called_from_GNRProd = true, est_df = colombianR_data, id = :id, time = :year, fixed_inputs = [:L :K], flexible_input = :RI, fes_returns = fes_res, int_const_series_degree = 3, lm_tfp_degree = 3, opts = opts);


gnr_res = GNRProdEst.gnrprodest!(data = colombianR_data, 
            output = :RGO, 
            flexible_input = :RI, 
            fixed_inputs = [:L :K], 
            ln_share_flex_y_var = :share, 
            id = :id, 
            time = :year,
            opts = opts);


gnr_res = GNRProdEst.gnrprodest(data = colombianR_data, 
            output = :RGO, 
            flexible_input = :RI, 
            fixed_inputs = [:L :K], 
            ln_share_flex_y_var = :share, 
            id = :id, 
            time = :year,
            opts = opts);



# Test with Colombian data
colombian_data = CSV.read("C:/Users/marku/Documents/GNRProdEst/Other Programs/R-version/Columbia_311.csv", DataFrame)


gnrprod(output = "RGO", fixed = c("L", "K"), flex = "RI",
                   share = "share", id = "id", time = "year", data = data)


opts = Dict("fes_series_degree" => 2,
            "fes_method" => "OLS", # OLS is experimental. It gives (almost) the same results but the residual is slightly different. However, it does not seem to matter too much when you check the second stage results
            "fes_print_starting_values" => true,
            "fes_print_results" => true,
            "int_const_series_degree" => 2,
            "lm_tfp_degree" => 1,
            "ses_optimizer" => NelderMead(),
            "ses_optimizer_options" =>  Optim.Options(
                                       f_tol = 1e-9,
                                       x_tol = 1e-2,
                                       g_tol = 1e-10,),
            "ses_print_starting_values" => true,
            "ses_print_results" => true)


opts = Dict("fes_series_degree" => 2,
            "fes_method" => "OLS", # OLS is experimental. It gives (almost) the same results but the residual is slightly different. However, it does not seem to matter too much when you check the second stage results
            "fes_print_starting_values" => true,
            "fes_print_results" => true)

colombian_data.k = log.(colombian_data.K)
colombian_data.l = log.(colombian_data.L)
colombian_data.ri = log.(colombian_data.ri)


gnr_res = GNRProdEst.gnrprodest!(data = colombian_data, 
                                output = :RGO, 
                                flexible_input = :RI, 
                                fixed_inputs = [:K :L], 
                                ln_share_flex_y_var = :share, 
                                id = :id, 
                                time = :year,
                                opts = opts);
