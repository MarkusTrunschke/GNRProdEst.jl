## General estimation command
function gnrprodest(;data::DataFrame, 
                  output::Symbol, 
                  flexible_input::Symbol, 
                  fixed_inputs::Union{Symbol,Array{Symbol}}, 
                  ln_share_flex_y::Symbol = :NotDefinedByUser, 
                  id::Symbol, 
                  time::Symbol, 
                  fes_starting_values::Vector = [missing],
                  ses_starting_values::Vector = [missing],
                  share_degree::Int = 3,
                  lm_tfp_degree::Int = 3,
                  int_const_series_degree::Int = 3,
                  opts::Dict = Dict())

    est_data = copy(data)
    
    fes_returns, ses_returns = gnrprodest!(data = est_data, 
                                                 output = output, 
                                                 flexible_input = flexible_input, 
                                                 fixed_inputs = fixed_inputs, 
                                                 ln_share_flex_y = ln_share_flex_y, 
                                                 id = id, 
                                                 time = time, 
                                                 fes_starting_values = fes_starting_values,
                                                 ses_starting_values = ses_starting_values,
                                                 share_degree = share_degree,
                                                 lm_tfp_degree = lm_tfp_degree,
                                                 int_const_series_degree = int_const_series_degree,
                                                 opts = opts)
    
    return fes_returns, ses_returns, est_data
end

## General estimation command that modifies its inputs (mostly the dataframe)
function gnrprodest!(;data::DataFrame, 
                  output::Symbol, 
                  flexible_input::Symbol, 
                  fixed_inputs::Union{Symbol,Array{Symbol}}, 
                  ln_share_flex_y::Symbol = :NotDefinedByUser, 
                  id::Symbol, 
                  time::Symbol,
                  fes_starting_values::Vector = [missing],
                  ses_starting_values::Vector = [missing],
                  share_degree::Int = 3,
                  lm_tfp_degree::Int = 3,
                  int_const_series_degree::Int = 3,
                  opts::Dict = Dict()
                  )
              
    # Clean inputs to be of the correct types
    fixed_inputs, _ , fes_starting_values, ses_starting_values = GNR_input_cleaner!(fixed_inputs = fixed_inputs, flexible_input = flexible_input, fes_starting_values = fes_starting_values, ses_starting_values = ses_starting_values)

    # Throw error if user messed up
    error_throw_fnc_first_stage(data, output, flexible_input, fixed_inputs, ln_share_flex_y, opts)

    # Get additional options right
    opts = opts_filler(opts)

    ## Run data preparation to ensure inputs are working with internals
    _, _, ln_share_flex_y = prep_data!(data, output = output, flexible_input = flexible_input, fixed_inputs = fixed_inputs, ln_share_flex_y = ln_share_flex_y, id = id, time = time)

    ## Run first stage estimation
    fes_returns = gnrfirststage!(data = data, output = output, flexible_input = flexible_input, fixed_inputs = fixed_inputs, ln_share_flex_y = ln_share_flex_y, share_degree = share_degree, starting_values = fes_starting_values, opts = opts)

    ## Run second stage estimation
    ses_returns = gnrsecondstage!(data = data, id = id, time = time, fixed_inputs = fixed_inputs, flexible_input = flexible_input, starting_values = ses_starting_values, lm_tfp_degree = lm_tfp_degree, int_const_series_degree = int_const_series_degree, called_from_GNRProd = true, fes_returns = fes_returns, opts = opts)

    # Return return dictionarys from both stages
    return fes_returns, ses_returns
end

## First stage function
function gnrfirststage!(;data::DataFrame, output::Symbol, flexible_input::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, ln_share_flex_y::Symbol, share_degree::Int = 3, starting_values::Vector = [missing], opts::Dict=Dict())
    
    # Clean inputs to make them comform with the rest of the program
    fixed_inputs, flexible_input, starting_values, _ = GNR_input_cleaner!(fixed_inputs = fixed_inputs, flexible_input = flexible_input, fes_starting_values = starting_values)
 
    all_input_symbols::Vector{Symbol} = vec(hcat(fixed_inputs..., flexible_input))

    # Throw error if user messed up
    error_throw_fnc_first_stage(data, output, flexible_input, fixed_inputs, ln_share_flex_y, opts)

    # Get additional options right
    opts = opts_filler(opts)

    # Add constant to data if not yet in dataset
    if "constant" ∉ names(data)
        data.constant = ones(size(data)[1])
    end

    # Get Polynomial series
    polynom_input_symbols = polynom_series!(data = data, var_names = all_input_symbols, degree = share_degree)

    # Run first stage estimation
    γ_dash, fes_opt_res = fes_est(data = data, ln_share_flex_y = ln_share_flex_y, input_var_symbols = polynom_input_symbols, method = opts["fes_method"], starting_values = starting_values, opts = opts)
    
    # Calculate some quantities
    fes_res, E = fes_predictions!(data = data, ln_share_flex_y = ln_share_flex_y, flexible_input = flexible_input, input_var_symbols = polynom_input_symbols, γ_dash = γ_dash, output = output)
    
    # Print results
    if opts["fes_print_results"] == true
        fes_print_res(fes_res, opts)
    end

    # Put return objects in dictionary
    return_elements = Dict("γ" => fes_res.γ,
                           "γ_dash" => γ_dash,
                           "γ_flex" => fes_res.γ_flex,
                           "E" => E,
                           "polynom_series" => polynom_input_symbols,
                           "all_inputs" => all_input_symbols,
                           "fixed_inputs" => fixed_inputs,
                           "flexible_input" => flexible_input,
                           "share_degree" => share_degree,
                           "fes_optim_estimates" => fes_res,
                           "fes_optim_results" => fes_opt_res
                           )

    return return_elements
end

## First stage estimation function that does not modify inputs (wrapper that copies inputs before passing it on)
function gnrfirststage(;data::DataFrame, output::Symbol, flexible_input::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, ln_share_flex_y::Symbol, share_degree::Int = 3, starting_values::Vector = [missing], opts::Dict=Dict())

    # Copy data s.t. the program does not modify existing data
    data = copy(data)

    # Run input-modifying version
    return_elements = gnrfirststage!(data = data,
                                     output = output, 
                                     flexible_input = flexible_input,
                                     fixed_inputs = fixed_inputs, 
                                     ln_share_flex_y = ln_share_flex_y,
                                     share_degree = share_degree,
                                     starting_values = starting_values, 
                                     opts = opts)

    # Return first stage results and dataframe with new variables
    return return_elements, data
end

## First stage estimation function
function fes_est(; data::DataFrame, ln_share_flex_y::Symbol, input_var_symbols::Array{Symbol}, method::String, starting_values::Vector, opts::Dict)
    
    if method == "OLS" # The GNR replication code uses a non-linear regression of the shares onto the natural logarithm of the polynomials for the first stage. However, this seems unnecessary. I cannot see any reason not to simply take the exponential of the shares onto the polynomials. This is a linear regression estimated by OLS. It is much faster to calculate and much more robust because it is not a numerical optimization but has a simple analytical solution.
        # Define matrices for OLS
        Y = exp.(data[:,ln_share_flex_y])
        X = hcat(data.constant, Matrix(data[:,input_var_symbols]))

        # Calculate OLS
        fes_results = vec(inv(X'*X)*X'*Y)

        fes_res = [] # No optimizer results to return b/c analytical solution

    elseif method == "NLLS"

        # Get starting values for NLLS (Follow replication code of GNR for now)
        γ_start = startvalues(data = data, Y_var = ln_share_flex_y, X_vars = input_var_symbols, user_start_vals = starting_values, stage = "first stage", opts = opts)
        
        # Get matrices
        X = Matrix(hcat(data.constant, Matrix(data[:,input_var_symbols])))
        Y = data[:, ln_share_flex_y]
        # Y = exp.(data[:, ln_share_flex_y])

        if any((X*γ_start) .< 0)
            throw("Error: Innver formula of NLLS function is already negative at starting values. Provide starting values manually or check for errors in data.")
        end

        # Allocate cache
        c = (cache1 = Array{Float64}(undef, length(Y), 1),
                obj = [1.0]
        )
        
        fes_res = optimize(γ -> NLLS_criterion!(γ, Y, X, c), γ_start, opts["fes_optimizer"], opts["fes_optimizer_options"])

        # Stop if algorithm did not converge
        if !Optim.converged(fes_res)
            throw("First stage optimization algorithm did not converge! Choose another optimizer, change optimizer settings or convergence criteria, and inspect your data carefully.")
        end

        fes_results = fes_res.minimizer # Return parameters

    end

    # Return solution vector
    return fes_results, fes_res
end

## Nonlinear least squares criterion function for the first stage
function NLLS_criterion!(γ, Y, X, c)
    
    mul!(c.cache1,X,γ)
    c.cache1 .= (Y .- NaNMath.log.(c.cache1)) .^2

    return sum!(c.obj, c.cache1)[1]
end

## Function to calculate quantities in the first stage (after estimation)
function fes_predictions!(;data::DataFrame, ln_share_flex_y::Symbol, flexible_input::Symbol, input_var_symbols::Array{Symbol}, γ_dash::Vector{<:Number}, output)
    
    flex_elas_sym = Symbol(flexible_input, "_elas")
    
    # Define matrices for OLS
    X_symbols = vcat(:constant, input_var_symbols)
    X = Matrix(data[!, X_symbols]) # hcat(data.constant, Matrix(data[:,input_var_symbols]))

    ## Calculate flexible input elasticity
    data.ln_D_E = NaNMath.log.(X*γ_dash) # Eq. (21)
    data.D_E = X*γ_dash

    # Calculate first stage residual ϵ
    data.ϵ = data.ln_D_E .- data[:,ln_share_flex_y]  # GNR and R-version of GNR calculate the residual like Xβ-Y instaed of  Y-Xβ. This is because of equation (11). There it is share = D_E - ϵ but we estimate D_E + ϵ. So take the negative of ϵ
    
    # Calculate constant E
    E = mean(exp.(data.ϵ))

    # Correct γ estimates
    γ = γ_dash ./ E # Below Eq. (21)

    # Calculate elasticity of the flexible input
    data[!, flex_elas_sym] = data.D_E / E # Eq. (14) (same as X*γ)
    
    # Calculate the integral (This part does the same as the R command (gnrflex) but differs from GNR's replication code. They never correct their coefficients for E.)
    # Get degree of intermediate input
    flex_degree = get_input_degree(flexible_input, X_symbols)
    γ_flex = γ ./ ( 1 .+ flex_degree[2,:]) # Calculate first part of first/second equation on p.2995

    data.int_flex = X*γ_flex .* data[! ,flexible_input] # Multiply by flex input to add the + 1

    # Calculate \mathcal(Y) (From eq. (16) and hint on p.2995)
    data.mathcal_Y = data[!, output] - data.int_flex - data.ϵ

    # Put results in DataFrame
    fes_res_df = DataFrame(Variable = X_symbols, γ = γ, γ_dash = γ_dash, γ_flex = γ_flex)

    # Return results
    return fes_res_df, E
end

## Function to print first stage results
function fes_print_res(fes_res::DataFrame, opts::Dict)
    # Print results if user wants that
    print_tab = hcat(fes_res.Variable, fes_res.γ, fes_res.γ_dash)
    header = (["Variable", "γ", "γ'"])

    println("First stage results:")
    pretty_table(print_tab, header = header, formatters = ft_printf("%5.5f"), limit_printing = false)

    return nothing
end

## Second stage estimation function
function gnrsecondstage!(;data::DataFrame, flexible_input::Symbol, fixed_inputs::Union{Array{Symbol},Symbol}, id::Symbol, time::Symbol, fes_returns::Dict, mathcal_Y_var::Symbol = :mathcal_Y, int_const_series_degree::Int = 3, lm_tfp_degree::Int = 3, called_from_GNRProd::Bool = false, starting_values::Vector = [missing], opts::Dict= Dict())
    
    # Clean inputs to be of the correct types
    fixed_inputs, _ , _ , starting_values  = GNR_input_cleaner!(fixed_inputs = fixed_inputs, flexible_input = flexible_input, ses_starting_values = starting_values)

    # Throw error if user messed up
    error_throw_fnc_sec_stage(data, flexible_input, fixed_inputs, id, time, opts)

    # Get additional options right
    if called_from_GNRProd == false
        opts = opts_filler(opts)
    end

    # Add constant to data if not yet in dataset
    if "constant" ∉ names(data)
        data.constant = ones(size(data)[1])
    end

    # Get starting values for GMM (Follow replication code of GNR for now)
    fixed_poly = Array{Symbol}[]
    if called_from_GNRProd == true # If called from within GNRProd, polynomials already exist. Use pure fixed input polynomials (as in GNR replication code)
        if typeof(fixed_inputs) == Symbol && int_const_series_degree == 1 # If there is only one fixed input and the polynomail degree is one, there is no need to select the correct polynomials
            fixed_poly = fixed_inputs
        else
            s_vec = string.(fes_returns["polynom_series"])
            strings_to_check = string.(fixed_inputs)

            polynom_series_checked = check_array_string_only_substrings(s_vec = s_vec, strings_to_check = strings_to_check)
            fixed_poly = Symbol.(polynom_series_checked[findall(polynom_series_checked[:, 2]), :][:,1])
        end
    else # If not called from within GNRProd, generate polynomials from fixed inputs
        fixed_poly = polynom_series!(data = data, var_names = fixed_inputs, degree = int_const_series_degree)
    end
    
    ses_starting_values = startvalues(data = data, Y_var = :mathcal_Y, X_vars = fixed_poly, user_start_vals = starting_values, stage = "second stage", opts = opts)[2:end]
    
    # Calculate some lags
    panel_lag!(data = data, id = id, time = time, variable = [mathcal_Y_var, fixed_poly...], lag_prefix = "lag_", lags = 1, drop_missings = true, force = true)

    # Run estimation
    ses_gmm_res, δ = ses_est(;data = data, starting_values = ses_starting_values, fixed_poly = fixed_poly, mathcal_Y_var = mathcal_Y_var, lm_tfp_degree = lm_tfp_degree, opts = opts)

    # Calculate some quantities
    ses_predictions!(data = data, mathcal_Y_var = mathcal_Y_var, flexible_input = flexible_input, fixed_inputs = fixed_inputs, fixed_poly = fixed_poly, fes_returns = fes_returns, α = Optim.minimizer(ses_gmm_res))

    # Print results
    if opts["ses_print_results"] == true
        ses_print_res(data = data, all_inputs = vec(hcat(fixed_inputs..., flexible_input)), fixed_poly = fixed_poly, δ = δ, α = Optim.minimizer(ses_gmm_res), lm_tfp_degree =lm_tfp_degree, opts = opts)
    end

    # Put return objects into dictionary
    ses_returns = Dict("α" => Optim.minimizer(ses_gmm_res),
                       "δ" => δ,
                       "gmm_optim_results" => ses_gmm_res,
                       "gmm_optim_options" => opts["ses_optimizer_options"],
                       "polynom_fixed" => fixed_poly,
                       "fixed_inputs" => fixed_inputs,
                       "flexible_input" => flexible_input,
                       "int_const_series_degree" => int_const_series_degree,
                       "lm_tfp_degree" => lm_tfp_degree,
                       "all_inputs" => fes_returns["all_inputs"]
    )

    return ses_returns
end

## Second stage estimation function that does not modify inputs (wrapper that copies inputs before passing it on)
function gnrsecondstage(;data::DataFrame, flexible_input::Symbol, fixed_inputs::Union{Array{Symbol},Symbol}, id::Symbol, time::Symbol, fes_returns::Dict, mathcal_Y_var::Symbol = :mathcal_Y, int_const_series_degree::Int = 3, lm_tfp_degree::Int = 3, called_from_GNRProd::Bool = false, starting_values::Vector = [missing], opts::Dict= Dict())

   # Copy data s.t. the program does not modify existing data
   data_df = copy(data)

   # Run input-modifying version
   ses_quant_res = gnrsecondstage!(data = data_df,
                   flexible_input = flexible_input,
                   fixed_inputs = fixed_inputs,
                   id = id,
                   time = time,
                   fes_returns = fes_returns,
                   mathcal_Y_var = mathcal_Y_var,
                   int_const_series_degree = int_const_series_degree,
                   lm_tfp_degree = lm_tfp_degree, 
                   called_from_GNRProd = called_from_GNRProd, 
                   starting_values = starting_values, 
                   opts = opts)

    # Return second stage returns and modified dataframe
    return ses_quant_res, data_df
end

## Second stage estimation function
function ses_est(;data::DataFrame, fixed_poly::Vector{Symbol}, starting_values::Vector, mathcal_Y_var::Symbol, lm_tfp_degree::Int = 3, opts::Dict)

    ## Run GMM estimation
    # Preallocate inputs
    fixed_poly_mat = Array(data[!, fixed_poly])
    fixed_poly_lag_mat = Array(data[!, "lag_" .* string.(fixed_poly)])
    constant = vec(ones(size(fixed_poly_mat)[1],1))
    w_lag_mat = Array{Float64}(undef, length(constant), lm_tfp_degree)
    mathcal_Y_vec = vec(data[!, mathcal_Y_var])
    lag_mathcal_Y_vec = vec(data[!, Symbol("lag_",mathcal_Y_var)])

    # Prepare cache
    c = (X = hcat(constant, w_lag_mat),
         cache1 = Array{Float64}(undef, length(constant), 1),
         cache2 = Array{Float64}(undef, length(constant), 1),
         cache3 = Array{Float64}(undef, 1 + lm_tfp_degree, 1 + lm_tfp_degree),
         cache4 = Array{Float64}(undef, size(fixed_poly_lag_mat)),
         cache5 = Array{Float64}(undef, 1, size(fixed_poly_lag_mat)[2]),
         coefs  = Array{Float64}(undef, 1 + lm_tfp_degree, 1), # 1 + lm_tfp_degree b/c of constant
         criterion = Array{Float64}(undef, 1, 1)
    )

    # Minimize GMM criterion function
    gmm_res = optimize(par -> ses_gmm!(α = par, 
                                       mathcal_Y_vec = mathcal_Y_vec, 
                                       lag_mathcal_Y_vec = lag_mathcal_Y_vec, 
                                       fixed_poly_mat = fixed_poly_mat, 
                                       fixed_poly_lag_mat = fixed_poly_lag_mat, 
                                       lm_tfp_degree = lm_tfp_degree, 
                                       c = c), 
                       starting_values,
                       opts["ses_optimizer"], 
                       opts["ses_optimizer_options"])

    # Stop if algorithm did not converge
    if !Optim.converged(gmm_res)
        throw("Second stage optimization algorithm did not converge! Choose another optimizer, change optimizer settings or convergence criteria, and inspect your data carefully.")
    end

    # Calculate law of motion for Productivity
    # Calculate w and X
    mul!(c.cache1, fixed_poly_mat, Optim.minimizer(gmm_res))
    c.cache1 .= mathcal_Y_vec .- c.cache1 # c.cache1 is now w
    mul!(c.cache2, fixed_poly_lag_mat, Optim.minimizer(gmm_res))

    c.X[:,2] .= lag_mathcal_Y_vec .- c.cache2 # X = hcat(constant, w_lag_mat)   
    polynomial_fnc_fast!(@view(c.X[:,2:end]), lm_tfp_degree, par_cal = false) # Fill X with polynomials of w

    # OLS
    mul!(c.cache3, c.X', c.X) # X'X
    mul!(c.coefs, c.X', c.cache1) # X'*Y
    ldiv!(cholesky!(Hermitian(c.cache3)), c.coefs) # inv(X'*X)*X*y, see https://discourse.julialang.org/t/memory-allocation-left-division-operator/103111/3 for an explination why this is the fastest way (even though not optimal for ill-conditioned matrices)
    
    # Return GMM results
    return gmm_res, c.coefs
end

## GMM criterion function
function ses_gmm!(;α::Array{<:Number},
                   mathcal_Y_vec::Vector{<:Number},
                   lag_mathcal_Y_vec::Vector{<:Number}, 
                   fixed_poly_mat::Array{<:Number},
                   fixed_poly_lag_mat::Array{<:Number},
                   lm_tfp_degree::Int,
                   c::NamedTuple)

    # Calculate w and X
    mul!(c.cache1, fixed_poly_mat, α)
    c.cache1 .= mathcal_Y_vec .- c.cache1 # c.cache1 is now w
    mul!(c.cache2, fixed_poly_lag_mat, α)
    c.X[:,2] .= lag_mathcal_Y_vec .- c.cache2 # X = hcat(constant, w_lag_mat)

    polynomial_fnc_fast!(@view(c.X[:,2:end]), lm_tfp_degree, par_cal = false) # Fill X with polynomials of w

    # OLS
    mul!(c.cache3, c.X', c.X) # X'X
    mul!(c.coefs, c.X', c.cache1) # X'*Y
    ldiv!(cholesky!(Hermitian(c.cache3)), c.coefs) # inv(X'*X)*X*y, see https://discourse.julialang.org/t/memory-allocation-left-division-operator/103111/3 for an explination why this is the fastest way (even though not optimal for ill-conditioned matrices)

    # Calculate residuals
    mul!(c.cache2, c.X, c.coefs) # cache2 is now w_hat
    c.cache1 .= c.cache1 .- c.cache2 # c.cache1 is now residuals

    # Calculate moments
    c.cache4 .= fixed_poly_mat .* c.cache1
    sum!(c.cache5, c.cache4) # Vector of moments ( fixed_poly_mat .* c.cache1 =: Matrix of moment vectors where each column is one moment vector)
    c.cache5 .= c.cache5 ./ length(c.cache1)

    # Return
    return mul!(c.criterion,c.cache5, c.cache5')[1] # Return criterion value (moment_vector * transposed(moment_vector))
end

## Function to calculate some quantities with results from second stage GMM estimation
function ses_predictions!(;data::DataFrame, mathcal_Y_var::Symbol, flexible_input::Symbol, fixed_inputs::Union{Symbol,Vector{Symbol}}, fixed_poly::Vector{Symbol}, α::Vector{<:Number}, fes_returns::Dict )

    # log(TFP) and TFP
    data.ω = data[!, mathcal_Y_var] .- Array(data[!, fixed_poly])*α
    data.Ω = exp.(data.ω)
    data.v = data.ω .* data.ϵ # ϵ = residulas from first stage


    ## Fixed input elasticity
    # First part: ∂C(fixed_inputs)/∂fixed_inputs
    # 1. Match all polynomial series variables to one lower polynomial
    all_inputs = fes_returns["all_inputs"]
    series_degree = get_input_degree(all_inputs, fes_returns["polynom_series"]) # Degree of input in polynomial series

    # Preallocate output matrix
    constants = Array{Float64}(undef, length(data.ω), size(series_degree)[1] - 2)
    no_constants = Array{Float64}(undef, length(data.ω), size(series_degree)[1] - 2)

    j = 1 # Iterator over output matrix
    for der_var_row = 2:size(series_degree)[1] - 1 # eachrow(series_degree[begin + 1:end-1,:])

        dev_fix_series_degree = copy(series_degree) # Get one degree lower version of series
        dev_fix_series_degree[der_var_row,:] .= ifelse.(series_degree[der_var_row,:] .> 0, series_degree[der_var_row,:] .- 1, series_degree[der_var_row,:] ) # end-1 b/c only one flex input. If allows more inputs, change to end - number of flex inputs
        C_deg = dev_fix_series_degree[:,dev_fix_series_degree[end,:] .== 0] # Drop columns with zero in flex input row

        # Check now which columns of original poly series fits the derivative of the current row (variable)
        match_arr = Array{Union{Symbol, Int}}(undef, 2, size(C_deg)[2]) # This will be an array of original poly series symbols and the column index that fits in the derivative of the series
        match_arr[1,:] .= C_deg[1,:]
        i = 1
        for col in eachcol(C_deg[begin + 1:end,:]) # Leave first row out because it contains symbols of polynomials
            check_mat = series_degree[begin + 1:end,:] .== col # Check if elements of the current column fit elements of each column of the polynomial series
            check_vec = prod.(eachcol(check_mat)) # Resulting vector has 1 if all elements of a column fit the current checked column. Product of all elements of each column
            indices = findall(check_vec)# Get the indices (should only be one at a time) of the fitting column
            if size(indices) != (0,) # If one is found write it into output matrix
                match_arr[2,i] = indices[]
            else # If none is found return a zero
                match_arr[2,i] = 0
            end
            i += 1 # Increase column counter
        end
        
        # Select columns from data that correspond to the correct polynomials. If match_arr == 0, there needs to be a constant.
        deriv_C = Array{Float64}(undef, size(data[!, mathcal_Y_var])[1], size(match_arr)[2])
        for col in 1:size(match_arr,2)
            if match_arr[2,col] == 0
                deriv_C[:,col] .= 1
            else
                deriv_C[:,col] .= data[!,fes_returns["polynom_series"][match_arr[2,col]]]
            end
        end
        
        # Calculate value of part of the derivative of the polynomial series corresponding to the current variable (row in series_degree)
        constants[:, j] .= deriv_C * (series_degree[der_var_row,series_degree[end,:] .== 0] .* α)

        j += 1 # Increase loop counter
    end

    j = 1
    for der_var_row = 2:size(series_degree)[1] - 1 # eachrow(series_degree[begin + 1:end-1,:])

        dev_fix_series_degree = copy(series_degree) # Get one degree lower version of series
        dev_fix_series_degree[der_var_row,:] .= ifelse.(series_degree[der_var_row,:] .> 0, series_degree[der_var_row,:] .- 1, series_degree[der_var_row,:] ) # end-1 b/c only one flex input. If allows more inputs, change to end - number of flex inputs
        dev_fix_series_degree[end,:] .= dev_fix_series_degree[end,:] .+ 1
        
        noC_deg = copy(dev_fix_series_degree)

        # Check now which columns of original poly series fits the derivative of the current row (variable)
        match_arr = Array{Union{Symbol, Int}}(undef, 2, size(noC_deg)[2]) # This will be an array of original poly series symbols and the column index that fits in the derivative of the series
        match_arr[1,:] .= noC_deg[1,:]
        i = 1
        for col in eachcol(noC_deg[begin + 1:end,:]) # Leave first row out because it contains symbols of polynomials
            check_mat = series_degree[begin + 1:end,:] .== col # Check if elements of the current column fit elements of each column of the polynomial series
            check_vec = prod.(eachcol(check_mat)) # Resulting vector has 1 if all elements of a column fit the current checked column. Product of all elements of each column
            indices = findall(check_vec)# Get the indices (should only be one at a time) of the fitting column
            if size(indices) != (0,) # If one is found write it into output matrix
                match_arr[2,i] = indices[]
            else # If none is found return a zero
                match_arr[2,i] = 0
            end
            i += 1 # Increase column counter
        end

        # Select columns from data that correspond to the correct polynomials. If match_arr == 0, there needs to be a constant.
        deriv_noC = Array{Float64}(undef, size(data[!, mathcal_Y_var])[1], size(match_arr)[2])
        for col in 1:size(match_arr,2)
            if match_arr[2,col] == 0
                deriv_noC[:,col] .= 0
            else
                deriv_noC[:,col] .= data[!,fes_returns["polynom_series"][match_arr[2,col]]]
            end
        end

        # Calculate value of part of the derivative of the polynomial series corresponding to the current variable (row in series_degree)
        no_constants[:, j] .= deriv_noC * (series_degree[der_var_row, :] .* fes_returns["γ_flex"][begin + 1:end])

        
        j += 1 # Increase loop counter
    end

    # Calculate elasticities of fixed inputs (are the sum of the constant and non-constant elasticity parts)
    fixed_elas = no_constants .+ constants

    # Put results in dataframe
    j = 1
    for input in fixed_inputs
        data[!, Symbol(input, "_elas")] = fixed_elas[:,j] # Add elasticity vector with correct symbol to data
        j += 1 # Increase counter
    end

    # Return results
    return fixed_elas
end

## Function printing results of the second stage
function ses_print_res(; data::DataFrame, all_inputs::Array{Symbol}, fixed_poly, δ, α, lm_tfp_degree, opts::Dict)
    
    desc_table = Array{Union{Symbol,Float64}}(undef, length(all_inputs), 5)

    j = 1
    for inp in all_inputs
        desc_table[j,1] = inp
        desc_table[j,2] = mean(data[!, Symbol(inp,"_elas")])
        desc_table[j,3] = std(data[!, Symbol(inp,"_elas")])
        desc_table[j,4] = minimum(data[!, Symbol(inp,"_elas")])
        desc_table[j,5] = maximum(data[!, Symbol(inp,"_elas")])

        j += 1
    end

    prod_table = [:ω mean(data.ω) std(data.ω) minimum(data.ω) maximum(data.ω)
                  :Ω mean(data.Ω) std(data.Ω) minimum(data.Ω) maximum(data.Ω)
                  :v mean(data.v) std(data.v) minimum(data.v) maximum(data.v)]

    int_const_tab = [fixed_poly α]

    lm_prod_table = Array{Union{String, Char, <:Number}}(undef, lm_tfp_degree + 1, 2) # Because of constant
    for j = 1:lm_tfp_degree + 1
        if j > 2
            lm_prod_table[j,1] = "ω"*string(superscript_this!(string(j-1)))
        elseif j == 2
            lm_prod_table[j,1] = "ω"
        elseif j == 1
            lm_prod_table[j,1] = "constant"
        end

        lm_prod_table[j,2] = δ[j]
    end


    header_var = (["Variable", "Mean", "SD", "Min", "Max"])
    header_par = (["Variable", "Estimate"])

    # Print parameters of the GMM estimation
    println("Integration constant series parameters")
    pretty_table(int_const_tab, header = header_par, formatters =  ft_printf("%5.5f"), limit_printing = false)

    # Print summary stats of all output elasticities of inputs
    println("All output elasticities:")
    pretty_table(desc_table, header = header_var, formatters =  ft_printf("%5.5f"), limit_printing = false)

    # Print productivity stats and law of motion parameters
    println("Productivity:")
    pretty_table(prod_table, header = header_var, formatters =  ft_printf("%5.5f"), limit_printing = false)

    println("Productivity Law of Motion:")
    pretty_table(lm_prod_table, header = header_par, formatters =  ft_printf("%5.5f"), limit_printing = false)

    return nothing
end