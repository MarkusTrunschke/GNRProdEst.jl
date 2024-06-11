using GNRProdEst
using Test
using CSV
using DataFrames
using Statistics
using Profile, ProfileView
using Optim

# Test with R-version's Colombian test data (likely industry 311)
colombianR_data = CSV.read("/Users/markus_trunschke/Documents/GNRProdEst/Other Programs/R-version/colombian.csv", DataFrame)
colombianR_data.constant = ones(length(colombianR_data.RGO))

opts = Dict("fes_series_degree" => 3,
            "lm_tfp_degree" => 3,
            "int_const_series_degree" => 3,
            "fes_method" => "NLLS", # OLS is experimental. It gives (almost) the same results but the residual is slightly different. However, it does not seem to matter too much when you check the second stage results
            "fes_print_starting_values" => true,
            "fes_print_results" => true,
            "ses_print_starting_values" => true,
            "ses_print_results" => true)


# User inputs
fixed_inputs = [:K :L]
flexible_input = :RI
output = :RGO
id = :id
time = :year
data = colombianR_data
ln_share_flex_y_var = :share
share_degree = 3
lm_tfp_degree = 3
int_const_series_degree = 3

# gnrprodest default inputs
fes_starting_values::Vector = [missing]
ses_starting_values::Vector = [missing]


# Clean inputs to be of the correct types
fixed_inputs, flexible_input, all_inputs, fes_starting_values = GNRProdEst.GNR_input_cleaner!(fixed_inputs = fixed_inputs, flexible_input = flexible_input, all_inputs = Array{Symbol}(undef,1,1), fes_starting_values = fes_starting_values, ses_starting_values = ses_starting_values)

# Throw error if user messed up
GNRProdEst.error_throw_fnc(data, output, flexible_input, fixed_inputs, ln_share_flex_y_var, id, time, fes_starting_values, ses_starting_values, opts)

# Get additional options right
opts = GNRProdEst.opts_filler!(opts)

## Run data preparation to ensure inputs are working with internals
GNRProdEst.prep_data!(data, output = output, flexible_input = flexible_input, fixed_inputs = fixed_inputs, ln_share_flex_y_var = ln_share_flex_y_var, id = id, time = time)

## Run first stage estimation
fes_returns = GNRProdEst.gnrfirststage!(est_df = data, output = output, flexible_input = flexible_input, fixed_inputs = fixed_inputs, ln_share_flex_y_var = ln_share_flex_y_var, all_input_symbols = all_inputs, share_degree = share_degree, starting_values = fes_starting_values, opts = opts)

## Run second stage estimation
ses_returns = GNRProdEst.gnrsecondstage!(est_df = data, id = id, time = time, fixed_inputs = fixed_inputs, flexible_input = flexible_input, all_inputs = all_inputs, starting_values = ses_starting_values, lm_tfp_degree = lm_tfp_degree, int_const_series_degree = int_const_series_degree, called_from_GNRProd = true, fes_returns = fes_returns, opts = opts)
