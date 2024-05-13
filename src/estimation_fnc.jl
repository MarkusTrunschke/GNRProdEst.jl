## General estimation command
function GNRProd(;data::DataFrame, 
                  output::Symbol, 
                  flex_input::Symbol, 
                  fixed_inputs::Union{Symbol,Array{Symbol}}, 
                  ln_share_flex_y_var::Symbol = :NotDefinedByUser, 
                  id::Symbol, 
                  time::Symbol, 
                  fes_starting_values::Vector = [missing],
                  ses_starting_values::Vector = [missing],
                  opts::Dict = Dict())
    
    # Clean inputs to be of the correct types
    fixed_inputs = GNR_input_cleaner!(fixed_inputs = fixed_inputs, stage = 0)

    # Throw error if user messed up
    error_throw_fnc(data,output,flex_input,fixed_inputs,ln_share_flex_y_var,id,time,fes_starting_values,ses_starting_values,opts)

    # Get additional options right
    opts = opts_filler!(opts)

    ## Run data preparation to ensure inputs are working with internals
    est_df, all_var_symbols, all_input_symbols, ln_share_flex_y_var = prep_data!(data, output = output, flex_input = flex_input, fixed_inputs = fixed_inputs, ln_share_flex_y_var = ln_share_flex_y_var, id = id, time = time)

    ## Run first stage estimation
    fes_returns = GNRFirstStage(est_df = est_df, output = output, flex_input = flex_input, fixed_inputs = fixed_inputs, ln_share_flex_y_var = ln_share_flex_y_var, all_input_symbols = all_input_symbols, starting_values = fes_starting_values, opts = opts)

    ## Run second stage estimation
    ses_returns = GNRSecondStage(est_df = est_df, fes_returns = fes_returns, fixed_inputs = fixed_inputs, starting_values = ses_starting_values, called_from_GNRProd = true, opts = opts)

    return fes_returns, ses_returns, est_df
end
## First stage function
function GNRFirstStage(;est_df, output::Symbol, flex_input::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, ln_share_flex_y_var::Symbol, all_input_symbols::Array{Symbol}, starting_values::Vector = [missing], opts::Dict=Dict())
    
    # Get Taylor polynomials
    taylor_input_symbols = taylor_series!(data = est_df, var_names = all_input_symbols, order = opts["fes_series_order"])

    # Run first stage estimation
    γ_dash, fes_res_df = fes_est(data=est_df, ln_share_flex_y_var=ln_share_flex_y_var, input_var_symbols=taylor_input_symbols, method = opts["fes_method"], print_res = opts["fes_print_results"], starting_values = starting_values, opts = opts)
    
    # Calculate some quantities
    γ, γ_flex, E = fes_predictions!(data = est_df, ln_share_flex_y_var = ln_share_flex_y_var, flex_input = flex_input, input_var_symbols=taylor_input_symbols, γ_dash = γ_dash, output = output)

    # Put return objects in dictionary
    return_elements = Dict("γ" => γ,
                           "γ_dash" => γ_dash,
                           "γ_flex" => γ_flex,
                           "E" => E,
                           "taylor_series" => taylor_input_symbols,
                           "fes_estimates" => fes_res_df
                           )

    return return_elements
end

## First stage estimation function
function fes_est(; data::DataFrame, ln_share_flex_y_var::Symbol, input_var_symbols::Array{Symbol}, method::String, print_res::Bool, starting_values::Vector, opts::Dict)
    
    if method == "OLS" # The GNR replication code uses a non-linear regression of the shares onto the natural logarithm of the taylor polynomials for the first stage. However, this seems unnecessary. I cannot see any reason not to simply take the exponential of the shares onto the taylor polynomials. This is a linear regression estimated by OLS. It is much faster to calculate and much more robust because it is not a numerical optimization but has a simple analytical solution.
        # Define matrices for OLS
        Y = exp.(data[:,ln_share_flex_y_var])
        X = hcat(data.constant, Matrix(data[:,input_var_symbols]))

        # Calculate OLS
        fes_results = vec(inv(X'*X)*X'*Y)

    elseif method == "NLLS"

        # Get starting values for NLLS (Follow replication code of GNR for now)
        γ_start = startvalues(data = data, Y_var = ln_share_flex_y_var, X_vars = input_var_symbols, user_start_vals = starting_values, stage = "first stage", print_res = opts["fes_print_starting_values"])
    
        # Get matrices
        X = Matrix(hcat(data.constant, Matrix(data[:,input_var_symbols])))
        Y = data[:, ln_share_flex_y_var]

        # Define model
        function model(X, γ)
            return log.(X*γ)
        end

        # Solve for parameters
        fes_res = curve_fit(model, X, Y, γ_start)

        fes_results = fes_res.param

    end

    # Print results if user wants that
    all_symbols = vcat([:constant], input_var_symbols)
    if print_res == true
        print_tab = hcat(all_symbols, fes_results)
        header = (["Variable", "Value"])

        println("First stage results:")
        pretty_table(print_tab, header = header)
    end

    # Get results in DataFrame
    fes_res_df = DataFrame(Variable = all_symbols, Coefficient = fes_results)

    # Return solution vector
    return fes_results, fes_res_df
end

## Function to calculate quantities in the first stage (after estimation)
function fes_predictions!(;data::DataFrame, ln_share_flex_y_var::Symbol, flex_input::Symbol, input_var_symbols::Array{Symbol}, γ_dash::Vector{Float64}, output)
    
    # Define matrices for OLS
    X = hcat(data.constant, Matrix(data[:,input_var_symbols]))
    
    # Calculate flexible input elasticity
    data.ln_D_E = log.(X*γ_dash) # Eq. (21)
    data.D_E = X*γ_dash

    # Calculate first stage residual ϵ
    data.ϵ = data.ln_D_E .- data[:,ln_share_flex_y_var]  # GNR and R-version of GNR calculate the residual like Xβ-Y instaed of  Y-Xβ. This is because of equation (11). There it is share = D_E - ϵ but we estimate D_E + ϵ. So take the negative of ϵ
    
    # Calculate constant E
    E = mean(exp.(data.ϵ))

    # Correct γ estimates
    γ = γ_dash ./ E # Below Eq. (21)

    # Calculate elasticity of the flexible input
    data.flex_elas = data.D_E / E # Eq. (14) (same as X*γ)
    
    # Calculate the integral (This part does the same as the R command (gnrflex) but differs from GNR's replication code. They never correct their coefficients for E.)
    # Get degree of intermediate input
    flex_degree = get_flex_degree(flex_input, vcat(:constant, input_var_symbols))
    γ_flex = γ ./ flex_degree[2,:] # Calculate first part of first/second equation on p.2995
    data.int_flex = X*γ_flex .* data[! ,flex_input] # Multiply by flex input to add the + 1

    # Calculate \mathcal(Y) (From eq. (16) and hint on p.2995)
    data.weird_Y = data[!, output] - data.int_flex - data.ϵ

    # Return results
    return γ, γ_flex, E
end

# Second stage estimation function
function GNRSecondStage(;est_df::DataFrame, fixed_inputs::Union{Array{Symbol},Symbol}, id::Symbol, time::Symbol, called_from_GNRProd::Bool = false, starting_values::Vector = [missing], fes_returns::Dict = Dict(), opts::Dict= Dict())

    # Clean inputs to be of the correct types
    fixed_inputs = GNR_input_cleaner!(fixed_inputs = fixed_inputs, stage = 2)

    # Get additional options right
    if called_from_GNRProd == false
        opts = opts_filler!(opts)
    end

    # Get starting values for GMM (Follow replication code of GNR for now)
    taylor_fixed = Array{Symbol}[]
    if called_from_GNRProd == true # If called from within GNRProd, taylor polynomials already exist. Use pure fixed input polynomials (as in GNR replication code)
        if typeof(fixed_inputs) == Symbol && opts["ses_series_order"] == 1 # If there is only one fixed input and the taylor order is one, there is no need to select the correct taylor polynomials
            taylor_fixed = fixed_inputs
        else
            s_vec = string.(fes_returns["taylor_series"])
            strings_to_check = string.(fixed_inputs)

            taylor_series_checked = check_array_string_only_substrings(s_vec = s_vec, strings_to_check = strings_to_check)
            taylor_fixed = Symbol.(taylor_series_checked[findall(taylor_series_checked[:, 2]), :][:,1])
        end
    else # If not called from within GNRProd, generate taylor polynomials from fixed inputs
        taylor_fixed = taylor_series!(data = est_df, var_names = fixed_inputs, order = opts["ses_series_order"])
    end
    
    ses_starting_values = startvalues(data = est_df, Y_var = :weird_Y, X_vars = taylor_fixed, user_start_vals = starting_values, stage = "secod stage", print_res =  opts["ses_print_starting_values"])

    # Run estimation
    ses_results = ses_est(;data = est_df, starting_values = ses_starting_values, id = id, time = time, taylor_fixed = taylor_fixed)

    return ses_results
end

# # Second stage estimation function
function ses_est(;data::DataFrame, id::Symbol, time::Symbol, weird_Y_var::Symbol = :weird_Y, taylor_fixed::Vector{Symbol}, starting_values::Vector)
    
    # Calculate some lags
    panel_lag!(data = data, id = id, time = time, variable = [weird_Y_var, taylor_fixed...], lag_prefix = "lag_", lags = 1, drop_missings = false, force = true)

    # Run GMM estimation
    gmm_res = ses_gmm(data, "lag_" .* string.(taylor_fixed))

    return gmm_res
end

# 2nd stage GMM function
function ses_gmm(data::DataFrame, α::Union{Vector, Symbol}, lag_taylor_fixed, )
    w = data.weird_Y .- data[!, taylor_fixed] * α  # - ak*k - ak2*k2
    w_1  = data.lag_weird_Y - data[!, lag_taylor_fixed] * α # ak*k_1-ak2*k2_1
	w2_1 = w_1^2 # w_1*w_1 
	w3_1 = w_1^3 # w2_1*w_1

    # OLS
    X = hcat(w_1, w2_1, w3_1)
    coefs = vec(inv(X'*X)*X'*w)
    residuals = w - x*coefs

    # Calculate moments
    m_vec = data[!, taylor_fixed] * residuals
    # m1 = k*residuals
    # m2 = k2*residuals

    return w
end