## General estimation command
function GNRProd(;data::DataFrame, output::Symbol, flex_inputs::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, opts::Dict)
    
    ## Run data preparation to ensure inputs are working with internals
    est_df, all_var_symbols, all_input_symbols = prep_data(data, output = output, flex_inputs = flex_inputs, fixed_inputs = fixed_inputs)

    ## Run first stage estimation
    sol = GNRFirstStage(;est_df, output = output, flex_inputs = flex_inputs, fixed_inputs = fixed_inputs, share = :share_m_y, all_input_symbols = all_input_symbols, opts = opts)


    return sol
end
## First stage function
function GNRFirstStage(;est_df, output::Symbol, flex_inputs::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, share::Symbol, all_input_symbols::Array{Symbol}, opts::Dict=Dict("fs_series_order" => 2))

    # Get Taylor polynomials
    taylor_input_symbols = taylor_series!(data = est_df, var_names = all_input_symbols, order = opts["fs_series_order"])

    # Get starting values for NLLS (Follow replication code of GNR for now)
    γ_start = fs_startvalues(data = est_df, share_var = share, input_var_symbols = taylor_input_symbols)

    # Run first stage estimation
    fs_res = fs_est(data=est_df, share_var=share, input_var_symbols=taylor_input_symbols, starting_values=γ_start,  method = opts["fs_method"])

    # ## Print results
    # GNRPrintres_fs()

    return fs_res
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
function fs_startvalues(;data::DataFrame, share_var::Symbol, input_var_symbols::Array{Symbol})

    # Define matrices for OLS
    Y = data[:,share_var]
    X = hcat(data.constant, Matrix(data[:,input_var_symbols]))

    # Calculate OLS
    startvals = vec(inv(X'*X)*X'*Y)

    # Return result
    return startvals
end

## First stage estimation function
function fs_est(; data::DataFrame, share_var::Symbol, input_var_symbols::Array{Symbol}, starting_values::Array, method::String = NLLS)
    if method == "OLS" # The GNR replication code uses a non-linear regression of the shares onto the natural logarithm of the taylor polynomials for the first stage. However, this seems unnecessary. I cannot see any reason not to simply take the exponential of the shares onto the taylor polynomials. This is a linear regression estimated by OLS. It is much faster to calculate and much more robust because it is not a numerical optimization but has a simple analytical solution.
        # Define matrices for OLS
        Y = exp.(data[:,share_var])
        X = hcat(data.constant, Matrix(data[:,input_var_symbols]))

        # Calculate OLS
        fs_res = vec(inv(X'*X)*X'*Y)

    elseif method == "NLLS"

        @. model(x, p) = p[1] * x

        fs_res = curve_fit(model, x, y, [1.0])


    end

    # Return solution vector
    return fs_res
end


function start_val_calc(;data = [], revenue = "", free_inputs = [], fixed_inputs = [], startvals = [], prod_fnc_constant = prod_fnc_constant, printvals = false)

    all_inputs = vcat("l_".*free_inputs, fixed_inputs)
    if prod_fnc_constant == true
        all_inputs = vcat(all_inputs, "constant")
    end

    ninputs = size(all_inputs,1) # Number of inputs (including constant for prod. fnc. if specified)

    if startvals != [] # If starting values are given already by user, just hand them back and exit function

        # Check number of input variables
        if typeof(startvals) == DataFrame # If starting values are given as a dataframe (as in bootstrapping function), take entries corresponding to prod. fnc. variables from coefs column
            startvals = startvals[1:length(all_inputs), :Coefs]
        end
        startvals = startvals
        if printvals == true
            println("User provided starting values are "*string(startvals))
        end
    else # If no starting values are given by the user, run OLS and use results as starting values
        est_data = data[!,vcat(all_inputs, revenue)]
        est_data = dropmissing(est_data) # Drop missings
        input_vars = Matrix(est_data[!,all_inputs])
        
        startvals = vec(inv(input_vars'*input_vars)*input_vars'*Matrix(select(est_data,revenue)))
        if printvals == true
            println("Computed starting values are "*string(startvals))
        end
    end
    return startvals
end