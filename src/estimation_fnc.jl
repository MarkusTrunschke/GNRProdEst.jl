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
function fes_predictions!(;data::DataFrame, ln_share_flex_y_var::Symbol, flex_input::Symbol, input_var_symbols::Array{Symbol}, γ_dash::Vector{<:Number}, output)
    
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
    flex_degree = get_input_degree(flex_input, vcat(:constant, input_var_symbols))
    γ_flex = γ ./ ( 1 .+ flex_degree[2,:]) # Calculate first part of first/second equation on p.2995

    data.int_flex = X*γ_flex .* data[! ,flex_input] # Multiply by flex input to add the + 1

    # Calculate \mathcal(Y) (From eq. (16) and hint on p.2995)
    data.weird_Y = data[!, output] - data.int_flex - data.ϵ

    # Return results
    return γ, γ_flex, E
end

# Second stage estimation function
function GNRSecondStage(;est_df::DataFrame, flex_input::Symbol, fixed_inputs::Union{Array{Symbol},Symbol}, id::Symbol, time::Symbol, weird_Y_var::Symbol = :weird_Y, w_degree::Int = 3, called_from_GNRProd::Bool = false, starting_values::Vector = [missing], fes_returns::Dict = Dict(), opts::Dict= Dict())

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
    
    ses_starting_values = startvalues(data = est_df, Y_var = :weird_Y, X_vars = taylor_fixed, user_start_vals = starting_values, stage = "secod stage", print_res =  opts["ses_print_starting_values"])[2:end]
    
    # Calculate some lags
    est_df = panel_lag(data = est_df, id = id, time = time, variable = [weird_Y_var, taylor_fixed...], lag_prefix = "lag_", lags = 1, drop_missings = true, force = true)

    # Run estimation
    ses_gmm_res = ses_est(;data = est_df, starting_values = ses_starting_values, taylor_fixed = taylor_fixed, weird_Y_var = weird_Y_var, w_degree = w_degree, opts = opts)

    # Calculate some quantities
    α = [0.76657932, -0.05454093]
    ses_quant_res = ses_predictions!(data = est_df, weird_Y_var = weird_Y_var, flex_input = flex_input, fixed_inputs = fixed_inputs, taylor_fixed = taylor_fixed, α = α)

    return ses_quant_res
end

# Second stage estimation function
function ses_est(;data::DataFrame, taylor_fixed::Vector{Symbol}, starting_values::Vector, weird_Y_var::Symbol, w_degree::Int = 3, opts::Dict)

    ## Run GMM estimation
    # Preallocate inputs
    taylor_fixed_mat = Array(data[!, taylor_fixed])
    taylor_fixed_lag_mat = Array(data[!, "lag_" .* string.(taylor_fixed)])
    constant = vec(ones(size(taylor_fixed_mat)[1],1))
    w_lag_mat = Array{Float64}(undef, length(constant), w_degree)
    weird_Y_vec = vec(data[!, weird_Y_var])
    lag_weird_Y_vec = vec(data[!, Symbol("lag_",weird_Y_var)])

    # Prepare cache
    c = (X = hcat(constant, w_lag_mat),
         cache1 = Array{Float64}(undef, length(constant), 1),
         cache2 = Array{Float64}(undef, length(constant), 1),
         cache3 = Array{Float64}(undef, 1 + w_degree, 1 + w_degree),
         cache4 = Array{Float64}(undef, size(taylor_fixed_lag_mat)),
         cache5 = Array{Float64}(undef, 1, size(taylor_fixed_lag_mat)[2]),
         coefs  = Array{Float64}(undef, 1 + w_degree, 1), # 1+ w_degree b/c of constant
         criterion = Array{Float64}(undef, 1, 1)
    )

    # Minimize GMM criterion function
    gmm_res = optimize(par -> ses_gmm!(α = par, 
                                       weird_Y_vec = weird_Y_vec, 
                                       lag_weird_Y_vec = lag_weird_Y_vec, 
                                       taylor_fixed_mat = taylor_fixed_mat, 
                                       taylor_fixed_lag_mat = taylor_fixed_lag_mat, 
                                       w_degree = w_degree, 
                                       c = c), 
                       starting_values,
                       opts["ses_optimizer"], 
                       opts["ses_optimizer_options"])

    # Return GMM results
    return gmm_res
end

# GMM criterion function
function ses_gmm!(;α::Array{<:Number},
                   weird_Y_vec::Vector{<:Number},
                   lag_weird_Y_vec::Vector{<:Number}, 
                   taylor_fixed_mat::Array{<:Number},
                   taylor_fixed_lag_mat::Array{<:Number},
                   w_degree::Int,
                   c::NamedTuple)

    # Calculate w and X
    mul!(c.cache1, taylor_fixed_mat, α)
    c.cache1 .= weird_Y_vec .- c.cache1 # c.cache1 is now w
    mul!(c.cache2, taylor_fixed_lag_mat, α)
    c.X[:,2] .= lag_weird_Y_vec .- c.cache2 # X = hcat(constant, w_lag_mat)

    polynomial_fnc_fast!(@view(c.X[:,2:end]), w_degree, par_cal = false) # Fill X with polynomials of w

    # OLS
    mul!(c.cache3, c.X', c.X) # X'X
    mul!(c.coefs, c.X', c.cache1) # X'*Y
    ldiv!(cholesky!(Hermitian(c.cache3)), c.coefs) # inv(X'*X)*X*y, see https://discourse.julialang.org/t/memory-allocation-left-division-operator/103111/3 for an explination why this is the fastest way (even though not optimal for ill-conditioned matrices)

    # Calculate residuals
    mul!(c.cache2, c.X, c.coefs) # cache2 is now w_hat
    c.cache1 .= c.cache1 .- c.cache2 # c.cache1 is now residuals

    # Calculate moments
    c.cache4 .= taylor_fixed_mat .* c.cache1
    sum!(c.cache5, c.cache4) # Vector of moments ( taylor_fixed_mat .* c.cache1 =: Matrix of moment vectors where each column is one moment vector)
    c.cache5 .= c.cache5 ./ length(c.cache1)

    # Return
    return mul!(c.criterion,c.cache5, c.cache5')[1] # Return criterion value (moment_vector * transposed(moment_vector))
end

# Function to calculate some quantities with results from second stage GMM estimation
function ses_predictions!(;data::DataFrame, weird_Y_var::Symbol, flex_input::Symbol, fixed_inputs::Union{Symbol,Vector{Symbol}}, taylor_fixed::Vector{Symbol}, α::Vector{<:Number})
    
    # log(TFP) and TFP
    data.ω = data[!, weird_Y_var] .- Array(data[!, taylor_fixed])*α
    data.Ω = exp.(data.ω)
    data.v = data.ω .* data.ϵ # ϵ = residulas from first stage

    ## Fixed input elasticity
    # First part: ∂C(fixed_inputs)/∂fixed_inputs

    CONTINUE HERE: CHECK WHICH INPUT DEGREE REALLY NEEDED (first stage? includes flexible and fixed inputs in this part of R code?)
    all_inputs = hcat(fixed_inputs, flex_input)
    taylor_input_symbols = taylor_series!(data = est_df, var_names = all_inputs, order = opts["fes_series_order"])
    # print(taylor_input_symbols)
    # print(get_input_degree(fixed_input, vcat(taylor_input_symbols)))


    # Integration constant
    # C_coef <- constant_gmm$par # GMM results
    # constants <- lapply(1:(nrow(input_degree) - 1), FUN = function(i) {
    #   new_in_deg <- input_degree
    #   new_in_deg[i, ] <- ifelse(new_in_deg[i, ] > 0,
    #                             new_in_deg[i, ] - 1,
    #                             new_in_deg[i, ])
      
    #   new_C_deg <- new_in_deg[, new_in_deg[nrow(input_degree), ] == 0]
    #   C_match <- apply(new_C_deg, MARGIN = 2, FUN = match_gnr, degree_vec =
    #                      input_degree)
      
    #   deriv_C <- all_input[, C_match]
    #   deriv_C[is.na(deriv_C)] <- 1
    #   C <- deriv_C %*%
    #     t(t(input_degree[i, input_degree[nrow(input_degree), ] == 0]) * C_coef)
    # })




    # elas_noC <- lapply(1:(nrow(orig_input_degree) - 1), FUN = function(i) {
    #     new_in_deg <- orig_input_degree
    #     new_in_deg[i, ] <- ifelse(new_in_deg[i, ] > 0,
    #                               new_in_deg[i, ] - 1,
    #                               new_in_deg[i, ])
    
    #     new_in_deg[nrow(new_in_deg), ] <- new_in_deg[nrow(new_in_deg), ] + 1
    
    #     in_match <- apply(new_in_deg, MARGIN = 2, FUN = match_gnr, degree_vec =
    #                         orig_input_degree)
    
    #     deriv_input <- orig_input[, in_match]
    #     deriv_input[is.na(deriv_input)] <- 0
    #     elas <- deriv_input %*% t(t(orig_input_degree[i, ]) * (object$arg$D_coef))
    #   })
    #   elas = lapply(1:length(elas_noC), FUN = function(x) {
    #     elas_noC[[x]] + constants[[x]]
    #   })

    return data
end