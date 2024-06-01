## Function that converts inputs into correct types for my program
function GNR_input_cleaner!(;fixed_inputs::Union{Array{Symbol},Symbol}, flex_input::Union{Array{Symbol},Symbol}, all_inputs::Array{Symbol}, stage::Int)
    # Convert fixed input into vector of symbols if an array was given
    if typeof(fixed_inputs) == Array || typeof(fixed_inputs) == Matrix{Symbol}
        fixed_inputs = vec(fixed_inputs)
    end

    # # Convert vector of symbols into symbol if there is only one element
    # if typeof(fixed_inputs) == Vector && size(fixed_inputs) == (1,)
    #     fixed_inputs = fixed_inputs[1]
    # end

    if isdefined(all_inputs,1) == false # Only check first element but that should be enough to see if user provided the input
        all_inputs = [x for sublist in [flex_input, fixed_inputs] for x in (sublist isa Vector ? sublist : [sublist])]
    end

    if stage == 0 # If it runs from the highest level command

    elseif stage == 1 # For inputs that are only present in the first stage

    elseif stage == 2 # For inputs that are only present in the second stage

    end
    
    return fixed_inputs, flex_input, all_inputs
end

## Function that checks if every input makes sense and thows an error if the user messed up
function error_throw_fnc(data::DataFrame, 
    output::Symbol, 
    flex_input::Symbol, 
    fixed_inputs::Union{Symbol,Array{Symbol}}, 
    ln_share_flex_y_var::Symbol, 
    id::Symbol, 
    time::Symbol, 
    fes_starting_values::Vector{<:Number},
    ses_starting_values::Vector{<:Number},
    opts::Dict)

    # Check if variables are actually in DataFrame
    missing_str = string()
    all_var_symbols = [x for sublist in [output, flex_input, fixed_inputs, id, time] for x in (sublist isa Vector ? sublist : [sublist])]

    for var in all_var_symbols
        if string(var) ∉ names(data)
            missing_str = missing_str*string(var)
        end
    end
    if missing_str != ""
        throw("The following columns are not in the specified dataframe: "*missing_str)
    end
    
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
    if "ses_series_order" ∉ keys(opts)
        opts["ses_series_order"] = 2
    end
    if "ses_optimizer" ∉ keys(opts)
        opts["ses_optimizer"] = NelderMead()
    end
    if "ses_optimizer_options" ∉ keys(opts)
        opts["ses_optimizer_options"] = Optim.Options(iterations = 100000,
                                                      f_tol = 1e-9,
                                                      x_tol = 1e-12,
                                                      g_tol = 1e-13, # √(Σ(yᵢ-ȳ)²)/n ≤ 1.0e-13 (only sets g_abstol, not outer_g_abstol)
                                                      allow_f_increases = true,
                                                      show_trace = false,
                                                      extended_trace = false,
                                                      show_every = 1,
                                                      time_limit = NaN,
                                                      store_trace = false)
    end
    if "ses_print_starting_values" ∉ keys(opts) 
        opts["ses_print_starting_values"] => false
    end
    if "ses_print_results" ∉ keys(opts)
        opts["ses_print_results"] => false
    end
    
    return opts
end

## Function that checks a string to only contain defined substrings
function check_str_only_def_substr(s::String, strings_to_check::Union{Vector,String,Char})
    parts = split(s, "_")

    check_vec = falses(length(parts))
    i = 1

    for part in parts
        for stri in strings_to_check
            if part == string(stri) # Need to convert stri in string to make sure that if stri is lenght of 1 it works. Otherwise it would compare a string to a Char and return false even if content is the same
                check_vec[i] = true
                break
            end
        end
        i += 1
    end

    return all(check_vec)
end

## Function that iterates over an array of strings and checks if it only contains defined substrings (or "_")
function check_array_string_only_substrings(;s_vec::Vector{String}, strings_to_check::Union{Vector,String,Char})
    check_res = falses(length(s_vec))
    j = 1
    for stri in s_vec
        res = check_str_only_def_substr(stri, [strings_to_check..., '_'])
        check_res[j] = res
        j += 1
    end
    return hcat(s_vec,check_res)
end

## Function to determine the intermediate input variable degree of the polynomial approximation based on the internal naming scheme
function get_input_degree(input::Union{Symbol,Vector{Symbol}}, all_var_symbols::Union{Symbol,Vector{Symbol}})
    
    # Clean inputs
    if typeof(input) == Symbol
        input = [input]
    end

    # Initalize array
    input_degree_mat = Array{Union{Symbol,<:Number}}(undef, length(input) + 1, length(all_var_symbols))# Array{Union{Symbol,Int}}[]
    
    # Give first row all variable symbols
    input_degree_mat[1,:] = all_var_symbols

    i = 2
    for inp in input # Loop over all inputs
        j = 1
        for sym in all_var_symbols # Iterate over all symbols of the Taylor polynomials
            parts = split(string(sym), "_")
            count_input = count(isequal(string(inp)), parts) # Check for each part if it matches the flexible input symbol and count the occurances
            input_degree_mat[i,j] = count_input

            j += 1 # Increase counter
        end
        i += 1 #  Increase counter
    end
    # Return matrix of symbols with count of input occurances
    return input_degree_mat
end

## Start values calculation function for the first stage (Currently only an OLS as in GNR 2020))
function startvalues(;data::DataFrame, Y_var::Symbol, X_vars::Union{Symbol,Array{Symbol}}, user_start_vals::Vector{<:Union{Missing, Number}}, stage::String, print_res::Bool)
    
    if any(ismissing.(user_start_vals)) || size(user_start_vals) != (1,) # If user did not specify starting values

        # Define matrices for OLS
        Y = data[:, Y_var]
        X = hcat(data.constant, Matrix(data[:,X_vars]))

        # Calculate OLS
        startvals = vec(inv(X'*X)*X'*Y)

        if stage == "first stage" # GNR make a correction to the constant in their first stage (this is probably because the constant is so strongly negative for their example that the ln(Xγ) returns NaNs because of negative Xγ-values)
            startvals[1] = 0.1 - minimum((X*startvals) .- startvals[1])
        end
    else # If user specified starting values
        
        # Check dimensions and throw an error if they do not match
        if size(X_vars) != size(user_start_vals)
            println(size(X_vars))
            println(size(user_start_vals))
            throw("You specified either too many or not enough "*stage*" starting values!")

        else
            startvals = user_start_vals
        end
    end

    # Print results if specified
    if print_res == true
        print_tab = hcat(vcat([:constant], X_vars), startvals)
        header = (["Variable", "Value"])
        println("First "*stage*" starting values:")
        pretty_table(print_tab, header = header)
    end

    # Return result
    return startvals
end

## Taylor Approximation function
function taylor_series!(;data::DataFrame, var_names::Union{Vector{Symbol},Symbol}, order::Int)

    if order > 4
        throw(error("Taylor series order is too large. The program only supports Taylor series until order 4"))
    end

    if typeof(var_names) == Symbol
        var_names = [var_names]
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

## Polynomial generation function
function polynomial_fnc!(data::DataFrame, variable::Symbol, degree::Int; name::Symbol = variable, force::Bool = false)
    
    # Start polynomial names list
    polynomials::Vector = [variable]

    # Calculate polynomials
    for i in 2:degree
        
        # Define polynomial column name
        new_sym = Symbol(name, "_", i)

        # Check if conflicting columns exist
        if force == false && string(new_sym) ∈ names(data)
            throw("Column name "*string(new_sym)*" already exists in dataframe. Either specify a different column symbol using the name option, rename conflicting columns in the dataframe, or use force option.")
        end

        # Calculate polynomial and put it in dataframe
        data[!, new_sym] = data[!, variable] .^ i

        # Add variable name to list
        push!(polynomials, new_sym)
    end

    return polynomials, data
    
end

## If there are already prepared columns in the dataframe and their values just need to be updated, jump in here. This is the fast version with no dynamic allocations at runtime.
function polynomial_fnc_fast!(poly_mat::Union{Array{<:Number}, SubArray{<:Number}}, degree::Int; par_cal::Bool = false)
    # Compute polynomial columns (each column of the matrix represents the i's polynomial of the first column)
    if par_cal == false
        for i in 2:degree
            poly_mat[:,i] .= @view(poly_mat[:,1]) .^ i
        end
    else
        Threads.@threads for i in 2:degree
            poly_mat[:,i] .= @view(poly_mat[:,1]) .^ i
        end
    end

    return poly_mat 
end

## Function that generates lagged values in a panel
# Panel lag function with return df
function panel_lag(;data::DataFrame, id::Symbol, time::Symbol, variable::Union{Symbol,Vector{Symbol}}, lag_prefix::String = "lag_", lags::Int = 1, drop_missings::Bool = false, force::Bool = false)
    
    # Clean input
    if typeof(variable) == Symbol
        variable = [variable]
    end

    # Sort data and create new df. This is the only difference to panel_lag!()
    new_df = copy(sort(data, [id, time]))

    # Drop lag columns if they already exist if force is true. Throw error if not.
    for var in variable
        if string(lag_prefix*string(var)) ∈ names(data)
            if force == true
                new_df = new_df[!, Not(Symbol(lag_prefix*string(var)))]
            else
                throw("Specified name for lag of variable already present in specified dataframe. Either set force = true, choose difference lag variable name, or rename the column.")
            end
        end
    end

    # Do the actual lagging
    lagging_that_panel!(data = new_df, id = id, time = time, variable = variable, lag_prefix = lag_prefix, lags = lags, drop_missings = drop_missings)
    
    # Return resulting Dataframe
    return new_df
end

# Panel lag function manipulating the original df
function panel_lag!(;data::DataFrame, id::Symbol, time::Symbol, variable::Union{Array{Symbol},Symbol}, lag_prefix::String = "lag_", lags::Int = 1, drop_missings::Bool = false, force::Bool = false)
    
    # Clean input
    if typeof(variable) == Symbol
        variable = [variable]
    end

    # Drop lag columns if they already exist if force is true. Throw error if not.
    for var in variable
        if string(lag_prefix*string(var)) ∈ names(data)
            if force == true
                data = data[!, Not(Symbol(lag_prefix*string(var)))]
            else
                throw("Specified name for lag of variable already present in specified dataframe. Either set force = true, choose difference lag variable name, or rename the column.")
            end
        end
    end

    # Sort data
    sort!(data, [id, time])

    # # Do the actual lagging
    data = lagging_that_panel!(data = data, id = id, time = time, variable = variable, lag_prefix = lag_prefix, lags = lags, drop_missings = drop_missings)

    return data
end

function lagging_that_panel!(;data::DataFrame, id::Symbol, time::Symbol, variable::Union{Symbol,Vector{Symbol}}, lag_prefix::String = "lag_", lags::Int = 1, drop_missings::Bool = false)

    # Generate lagged values per id and select all but the original variable (causes problems in join). The ShiftedArrays.lag function names the lagged column itself with variable_lag
    df_lag = select(data, [id, time, variable...])
    
    lag_variables::Vector{Symbol} = []

    # Lag all variables
    for lag_var in [time, variable...]
        transform!(groupby(df_lag, id), lag_var => ShiftedArrays.lag)

        push!(lag_variables, Symbol(string(lag_var) .* "_lag"))
    end

    # Drop missings in lagged variables we just generated
    dropmissing!(df_lag, lag_variables)
    
    # Check if lag is actually only the expected lag apart
    for var in lag_variables[2:end]
        df_lag[!, var] = ifelse.(df_lag[!,time] .- lags .== df_lag[!,Symbol(string(time)*"_lag")], df_lag[!,var], missing)
    end

    select!(df_lag, [time, id, lag_variables[2:end]...]) # Drop lagged time variable from df

    # Combine lagged variable with original data and sort it.
    sort!(leftjoin!(data, df_lag, on = [id, time]), [id, time])

    # Drop missings in lagged variable we just generated if user wants to
    if drop_missings == true
        dropmissing!(data, lag_variables[2:end])
    end
    
    # Rename variable to user-specified name
    for var in variable
        rename!(data, Symbol(string(var)*"_lag") => Symbol(lag_prefix*string(var)))
    end

    # Return result
    return data
end