## General estimation command
function GNRProd(;data::DataFrame, output::Symbol, flex_input::Symbol, fixed_inputs::Union{Symbol,Array{Symbol}}, ln_share_flex_y_var::Symbol = :NotDefinedByUser, opts::Dict = Dict())
    
    # Get additional options right
    opts = opts_filler!(opts)

    ## Run data preparation to ensure inputs are working with internals
    est_df, all_var_symbols, all_input_symbols, ln_share_flex_y_var = prep_data!(data, output = output, flex_input = flex_input, fixed_inputs = fixed_inputs, ln_share_flex_y_var = ln_share_flex_y_var)

    ## Run first stage estimation
    fes_returns = GNRFirstStage(est_df = est_df, output = output, flex_input = flex_input, fixed_inputs = fixed_inputs, ln_share_flex_y_var = ln_share_flex_y_var, all_input_symbols = all_input_symbols, opts = opts)

    ## Run second stage estimation
    ses_returns = GNRSecondStage()

    return fes_returns, ses_returns
end
## First stage function
function GNRFirstStage(;est_df, output::Symbol, flex_input::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, ln_share_flex_y_var::Symbol, all_input_symbols::Array{Symbol}, opts::Dict=Dict())
    
    # Get Taylor polynomials
    taylor_input_symbols = taylor_series!(data = est_df, var_names = all_input_symbols, order = opts["fes_series_order"])

    # Run first stage estimation
    γ_dash, fes_res_df = fes_est(data=est_df, ln_share_flex_y_var=ln_share_flex_y_var, input_var_symbols=taylor_input_symbols, method = opts["fes_method"], print_res = opts["fes_print_results"], opts = opts)
    
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
function fes_startvalues(;data::DataFrame, ln_share_flex_y_var::Symbol, input_var_symbols::Array{Symbol}, print_res::Bool)
    
    println(names(data))
    display(data[:, ln_share_flex_y_var])
    # Define matrices for OLS
    Y = data[:, ln_share_flex_y_var]
    X = hcat(data.constant, Matrix(data[:,input_var_symbols]))

    # Calculate OLS
    startvals = vec(inv(X'*X)*X'*Y)

    # GNR make a correction to the constant (this is probably because the constant is so strongly negative for their example that the ln(Xγ) returns NaNs because of negative Xγ-values)
    startvals[1] = 0.1 - minimum((X*startvals) .- startvals[1])

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
function fes_est(; data::DataFrame, ln_share_flex_y_var::Symbol, input_var_symbols::Array{Symbol}, method::String, print_res::Bool, opts::Dict)
    
    if method == "OLS" # The GNR replication code uses a non-linear regression of the shares onto the natural logarithm of the taylor polynomials for the first stage. However, this seems unnecessary. I cannot see any reason not to simply take the exponential of the shares onto the taylor polynomials. This is a linear regression estimated by OLS. It is much faster to calculate and much more robust because it is not a numerical optimization but has a simple analytical solution.
        # Define matrices for OLS
        Y = exp.(data[:,ln_share_flex_y_var])
        X = hcat(data.constant, Matrix(data[:,input_var_symbols]))

        # Calculate OLS
        fes_results = vec(inv(X'*X)*X'*Y)

    elseif method == "NLLS"

        # Get starting values for NLLS (Follow replication code of GNR for now)
        γ_start = fes_startvalues(data = data, ln_share_flex_y_var = ln_share_flex_y_var, input_var_symbols = input_var_symbols, print_res = opts["fes_print_starting_values"])

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
    ln_flex_share_y = data[:,ln_share_flex_y_var]
    X = hcat(data.constant, Matrix(data[:,input_var_symbols]))
    
    # Calculate flexible input elasticity
    # data.flex_elas = X*γ
    data.ln_D_E = X*γ_dash # GNR and R-version of GNR calculate this without correcting for the constant. Follow them

    # Calculate first stage residual ϵ
    data.ϵ = data.ln_D_E .- ln_flex_share_y  # GNR and R-version of GNR calculate the residual like Xβ-Y instaed of  Y-Xβ. This is because of equation (11). There it is share = D_E - ϵ but we estimate D_E + ϵ. So take the negative of ϵ
    
    # Calculate constant E
    E = mean(exp.(data.ϵ))

    γ = γ_dash ./ E # Below Eq. (21)

    data.flex_elas = exp.(data.ln_D_E) / E # Eq. (14)
    
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
function GNRSecondStage(;data::DataFrame)
    
    # Run estimation
    ses_results = ses_est(;data::DataFrame, )

end

# Second stage estimation function
function ses_est()
    
end