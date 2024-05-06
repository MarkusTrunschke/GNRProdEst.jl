## General estimation command
function GNRProd(;data::DataFrame, 
                  output::Symbol, 
                  flex_input::Symbol, 
                  fixed_inputs::Union{Symbol,Array{Symbol}}, 
                  ln_share_flex_y_var::Symbol = :NotDefinedByUser, 
                  id::Symbol, 
                  time::Symbol, 
                  fes_starting_values::Array{Union{Float,Int,Missing}} = [Missing],
                  ses_starting_values::Array{Union{Float,Int,Missing}} = [Missing],
                  opts::Dict = Dict())
    
    # Get additional options right
    opts = opts_filler!(opts)

    ## Run data preparation to ensure inputs are working with internals
    est_df, all_var_symbols, all_input_symbols, ln_share_flex_y_var = prep_data!(data, output = output, flex_input = flex_input, fixed_inputs = fixed_inputs, ln_share_flex_y_var = ln_share_flex_y_var, id = id, time = time)

    ## Run first stage estimation
    fes_returns = GNRFirstStage(est_df = est_df, output = output, flex_input = flex_input, fixed_inputs = fixed_inputs, ln_share_flex_y_var = ln_share_flex_y_var, all_input_symbols = all_input_symbols, starting_values = fes_starting_values, opts = opts)

    ## Run second stage estimation
    ses_returns = GNRSecondStage(est_df = est_df, fes_returns = fes_returns, called_from_GNRProd = true)

    return fes_returns, ses_returns
end
## First stage function
function GNRFirstStage(;est_df, output::Symbol, flex_input::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, ln_share_flex_y_var::Symbol, all_input_symbols::Array{Symbol}, starting_values::Array{Union{Float,Int,Missing}}, opts::Dict=Dict())
    
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

## Auxiliary function to fill up options that were not given in opts dictionary
function opts_filler!(opts::Dict)
    
    if "fes_series_order" ∉ keys(opts)
        opts["fes_series_order"] = 2
    end
    if "fes_print_starting_values" ∉ keys(opts)
        opts["fes_print_starting_values"] = false
    end
    if "fes_print_results" ∉ keys(opts)
        opts["fes_print_results"] = false
    end
    if "fes_method" ∉ keys(opts)
        opts["fes_method"] = "OLS"
    end
    if "ses_print_starting_values" ∉ keys(opts)
        opts["ses_print_starting_values"] = false
    end
    return opts
end

## Taylor Approximation function
function taylor_series!(;data::DataFrame, var_names::Vector{Symbol}, order::Int)

    if order > 4
        throw(error("Taylor series order is too large. The program only supports Taylor series until order 4"))
    end

    nvar = length(var_names)

    taylor_var_names::Array{Symbol} = []

    # Generate all variables
    for j = 1 : nvar
        push!(taylor_var_names, var_names[j]) # Add variable name to taylor variable symbol list

        # If order of taylor approximation is larger than 1
        if order > 1
            for i = j : nvar
                varn = Symbol(String(var_names[j])*"_"*String(var_names[i])) # Generate variable name of new variable
                push!(taylor_var_names, varn) # Add variable name to taylor variable symbol list
                if varn ∈ names(data) # Check if variable is already in dataset. If yes: Drop it
                    data = select(data, Not(varn))
                end
                data[:, varn] = data[:, var_names[j]] .* data[:,var_names[i]] # Add variable to dataset

                # If order of taylor approximation is larger than 2
                if order > 2
                    for u = i : nvar
                        varn = Symbol(String(var_names[j])*"_"*String(var_names[i])*"_"*String(var_names[u])) # generate variable name of new variable
                        push!(taylor_var_names, varn) # Add variable name to taylor variable symbol list
                        if varn ∈ names(data) # Check if variable is already in dataset. If yes: Drop it
                            data = select(data, Not(varn))
                        end
                        data[:,varn] = data[:,var_names[j]].*data[:,var_names[i]].*data[:,var_names[u]]

                        # If order of taylor approximation is larger than 3
                        if order > 3
                            for r = u : nvar
                                varn = Symbol(String(var_names[j])*"_"*String(var_names[i])*"_"*String(var_names[u])*"_"*String(var_names[r])) # generate variable name of new variable
                                push!(taylor_var_names, varn) # Add variable name to taylor variable symbol list
                                if varn ∈ names(data) # Check if variable is already in dataset. If yes: Drop it
                                    data = select(data, Not(varn))
                                end
                                data[:,varn] = data[:,var_names[j]].*data[:,var_names[i]].*data[:,var_names[u]].*data[:,var_names[r]]

                            end
                        end
                    end
                end
            end
        end
    end

    # Return result
    return taylor_var_names
end

## Start values calculation function for the first stage (Currently only an OLS as in GNR 2020))
function startvalues(;data::DataFrame, Y_var::Symbol, X_vars::Union{Symbol,Array{Symbol}}, user_start_vals::Array{Union{Float,Int,Missing}}, stage::String, print_res::Bool)
    
    if user_start_vals != [Missing] # If user specified starting values

        # Check dimensions and throw an error if they do not match
        if size(X_vars) != size(user_start_vals)
            throw("You specified either too many or not enough "*stage*" starting values!")

        else
            startvals = user_start_vals
        end

    else # If user did not specify starting values

        # Define matrices for OLS
        Y = data[:, Y_var]
        X = hcat(data.constant, Matrix(data[:,X_vars]))

        # Calculate OLS
        startvals = vec(inv(X'*X)*X'*Y)

        # GNR make a correction to the constant (this is probably because the constant is so strongly negative for their example that the ln(Xγ) returns NaNs because of negative Xγ-values)
        startvals[1] = 0.1 - minimum((X*startvals) .- startvals[1])
    end

    # Print results if specified
    if print_res == true
        print_tab = hcat(vcat([:constant], input_var_symbols), startvals)
        header = (["Variable", "Value"])
        println("First stage starting values:")
        pretty_table(print_tab, header = header)
    end

    # Return result
    return startvals
end

## First stage estimation function
function fes_est(; data::DataFrame, ln_share_flex_y_var::Symbol, input_var_symbols::Array{Symbol}, method::String, print_res::Bool, starting_values::Array{Union{Float,Int,Missing}}, opts::Dict)
    
    if method == "OLS" # The GNR replication code uses a non-linear regression of the shares onto the natural logarithm of the taylor polynomials for the first stage. However, this seems unnecessary. I cannot see any reason not to simply take the exponential of the shares onto the taylor polynomials. This is a linear regression estimated by OLS. It is much faster to calculate and much more robust because it is not a numerical optimization but has a simple analytical solution.
        # Define matrices for OLS
        Y = exp.(data[:,ln_share_flex_y_var])
        X = hcat(data.constant, Matrix(data[:,input_var_symbols]))

        # Calculate OLS
        fes_results = vec(inv(X'*X)*X'*Y)

    elseif method == "NLLS"

        # Get starting values for NLLS (Follow replication code of GNR for now)
        γ_start = fes_startvalues(data = data, ln_share_flex_y_var = ln_share_flex_y_var, input_var_symbols = input_var_symbols, user_start_vals = starting_values, stage = "first stage", print_res = opts["fes_print_starting_values"])

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

# Function to determine the intermediate input variable degree of the polynomial approximation based on the internal naming scheme
function get_flex_degree(flex_input::Symbol, all_var_symbols::Array{Symbol})
    
    # Initalize array
    flex_degree_list = Array{Union{Symbol,Int}}[]

    for sym in all_var_symbols # Iterate over all symbols of the Taylor polynomials
        parts = split(string(sym), "_")
        count_flex = count(isequal(string(flex_input)), parts) # Check for each part if it matches the flexible input symbol and count the occurances
        push!(flex_degree_list, [sym, count_flex + 1])
    end

    # Return list of symbols with count of flexible input occurances
    return hcat(flex_degree_list...)
end



# Second stage estimation function
function GNRSecondStage(;est_df::DataFrame, fixed_inputs::Union{Array{Symbol},Symbol}, called_from_GNRProd::Bool, fes_returns::Dict = Dict(), opts::Dict())
    
    # Get starting values for GMM (Follow replication code of GNR for now)
    taylor_fixed = Array{Symbol}[]
    if called_from_GNRProd == true # If called from within GNRProd, taylor polynomials already exist. Use pure fixed input polynomials (as in GNR replication code)
        if typeof(fixed_inputs) == Symbol && opts["fes_series_order"] == 1 # If there is only one fixed input and the taylor order is one, there is no need to select the correct taylor polynomials
            taylor_fixed = fixed_inputs
        else
            taylor_series_checked = check_array_string_only_substrings(s_vec = string.(fes_returns["taylor_series"]), strings_to_check = string.(fixed_inputs))
            taylor_fixed = Symbol.(taylor_series_checked[findall(res[:, 2]), :][:,1])
        end
    else # If not called from within GNRProd, generate taylor polynomials from fixed inputs
        taylor_fixed = taylor_series!(data = est_df, var_names = fixed_inputs, order = opts["ses_series_order"])
    end
    
    ses_starting_values = startvalues(data = data, Y_var = :weird_Y, X_vars = taylor_fixed, user_start_vals = ses_starting_values, stage = "secod stage", print_res =  opts["ses_print_starting_values"])

    # Run estimation
    ses_results = ses_est(;data::DataFrame, ses_starting_values = ses_starting_values)

end

# Second stage estimation function
function ses_est(;data::DataFrame)
    
end