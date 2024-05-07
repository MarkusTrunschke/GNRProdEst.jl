## Function that converts inputs into correct types for my program
function GNR_input_cleaner!(;fixed_inputs::Union{Array{Symbol},Symbol}, stage::Int)
    # Convert fixed input into vector of symbols if an array was given
    if typeof(fixed_inputs) == Array || typeof(fixed_inputs) == Matrix{Symbol}
        fixed_inputs = vec(fixed_inputs)
    end

    # Convert vector of symbols into symbol if there is only one element
    if typeof(fixed_inputs) == Vector && size(fixed_inputs) == (1,)
        fixed_inputs = fixed_inputs[1]
    end
    
    if stage == 0 # If it runs from the highest level command

    elseif stage == 1 # For inputs that are only present in the first stage

    elseif stage == 2 # For inputs that are only present in the second stage

    end
    
    return fixed_inputs
end

## Function that checks if every input makes sense and thows an error if the user messed up
function error_throw_fnc(data::DataFrame, 
    output::Symbol, 
    flex_input::Symbol, 
    fixed_inputs::Union{Symbol,Array{Symbol}}, 
    ln_share_flex_y_var::Symbol, 
    id::Symbol, 
    time::Symbol, 
    fes_starting_values::Vector,
    ses_starting_values::Vector,
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
    return opts
end

## Function that checks a string to only contain defined substrings
function check_str_only_def_substr(s::String, strings_to_check::Union{Vector{String},Vector{Char},String})
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
function check_array_string_only_substrings(;s_vec::Vector{String}, strings_to_check::Union{Vector{String},Vector{Char},String})
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

## Start values calculation function for the first stage (Currently only an OLS as in GNR 2020))
function startvalues(;data::DataFrame, Y_var::Symbol, X_vars::Union{Symbol,Array{Symbol}}, user_start_vals::Vector, stage::String, print_res::Bool)
    
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

## Function that generates lagged values in a panel
# Second panel lag function
function panel_lag(;data::DataFrame, id::Symbol, time::Symbol, variable::Symbol, lag_prefix::String = "l_", lags::Int = 1, drop_missings::Bool = false)

    # Sort data
    sdf = copy(sort(data, [id, time]))

    # Generate lagged values per id and select all but the original variable (causes problems in join). The ShiftedArrays.lag function names the lagged column itself with variable_lag
    sdf = select(transform(groupby(sdf, id), variable => ShiftedArrays.lag, time => ShiftedArrays.lag), [Symbol(string(variable)*"_lag"), time, Symbol(string(time)*"_lag"), id])

    # Drop missings in lagged time variable we just generated
    sdf = dropmissing(sdf, Symbol(string(time)*"_lag"))

    # Check if lag is actually only the expected lag apart
    sdf[!,Symbol(string(variable)*"_lag")] = ifelse.(sdf[!,time] .- lags .== sdf[!,Symbol(string(time)*"_lag")], sdf[!,Symbol(string(variable)*"_lag")], missing)



    sdf = select(sdf, Not(string(time)*"_lag")) # Drop lagged time variable from df

    # Combine lagged variable with original data and sort it.
    newdata = sort(leftjoin(data, sdf, on = [id, time]), [:id, :time])

    # Drop missings in lagged variable we just generated if user wants to
    if drop_missings == true
        newdata = dropmissing(newdata, Symbol(string(variable)*"_lag"))
    end
    
    # Rename variable to user-specified name
    rename!(newdata, Symbol(string(variable)*"_lag") => Symbol(lag_prefix*string(variable)))

    # Return result
    return newdata
end