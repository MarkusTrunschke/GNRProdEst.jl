## Function that converts inputs into correct types for my program
function GNR_input_cleaner!(;fixed_inputs::Union{Array{Symbol},Symbol}, flexible_input::Union{Array{Symbol},Symbol}, fes_starting_values::Union{Vector{<:Number},Vector{Missing}} = vec([missing]), ses_starting_values::Union{Vector{<:Number},Vector{Missing}} = vec([missing]))
    # Convert fixed input into vector of symbols if an array was given
    if typeof(fixed_inputs) == Array || typeof(fixed_inputs) == Matrix{Symbol}
        fixed_inputs = vec(fixed_inputs)
    end

    # Convert fixed input to a vector if user put in a symbol. (Makes working with it in the program easier) 
    if typeof(fixed_inputs) == Symbol
        fixed_inputs = [fixed_inputs]
    end
    
    # Convert Int to Float because optimisers don't like Int starting values
    if typeof(fes_starting_values) != Vector{Missing} && size(fes_starting_values) != (1,)
        if typeof(fes_starting_values) == Vector{Int64}
            fes_starting_values = float.(fes_starting_values)
        end
    end

    # Convert Int to Float because optimisers don't like Int starting values
    if typeof(ses_starting_values) != Vector{Missing} && size(ses_starting_values) != (1,)
        if typeof(ses_starting_values) == Vector{Int64}
            ses_starting_values = float.(ses_starting_values)
        end
    end

    return fixed_inputs, flexible_input, fes_starting_values, ses_starting_values #, opts_any # Return cleaned inputs
end

## Function that checks if every input makes sense and thows an error if the user messed up
function error_throw_fnc(data::DataFrame, 
                         output::Symbol, 
                         flexible_input::Symbol, 
                         fixed_inputs::Union{Symbol,Array{Symbol}}, 
                         ln_share_flex_y_var::Symbol,
                         id::Symbol, 
                         time::Symbol,
                         opts::Dict)

    # Check if variables are actually in DataFrame
    missing_str = string()
    non_num_str = string()

    all_var_symbols = Array{Symbol}(undef,0)
    if ln_share_flex_y_var != :NotDefinedByUser
        all_var_symbols = vec(hcat(fixed_inputs..., flexible_input, id, time, output, ln_share_flex_y_var))
    else
        all_var_symbols = vec(hcat(fixed_inputs..., flexible_input, id, time, output))
    end

    for var in all_var_symbols
        if string(var) ∉ names(data)
            missing_str = missing_str*string(var)*", "
        end
    end
    if missing_str != ""
        throw("The following columns are not in the specified dataframe: "*missing_str[1:end-2])
    end

    # Check if input, output, and ln_share_flex_y are numeric
    for var in setdiff(all_var_symbols, vec(hcat(id, time)))
        if !(eltype(data[!, var]) <:Union{Missing, Number})
            non_num_str = non_num_str*string(var)*", "
        end
    end
    if non_num_str != ""
        throw("The following columns should be numeric but are not: "*non_num_str[1:end-2])
    end
    
end

function error_throw_fnc_first_stage(data::DataFrame, 
    output::Symbol,
    flexible_input::Symbol, 
    fixed_inputs::Union{Symbol,Array{Symbol}}, 
    ln_share_flex_y_var::Symbol,
    opts::Dict)

    # Check if variables are actually in DataFrame
    missing_str = string()
    non_num_str = string()
    
    all_var_symbols = Array{Symbol}(undef,0)
    if ln_share_flex_y_var != :NotDefinedByUser
        all_var_symbols = vec(hcat(fixed_inputs..., flexible_input, output, ln_share_flex_y_var))
    else
        all_var_symbols = vec(hcat(fixed_inputs..., flexible_input, output))
    end

    for var in all_var_symbols
        if string(var) ∉ names(data)
            missing_str = missing_str*string(var)*", "
        end
    end
    if missing_str != ""
        throw("The following columns are not in the specified dataframe: "*missing_str[1:end-2])
    end

    # Check if input, output, and ln_share_flex_y are numeric
    for var in all_var_symbols
        if !(eltype(data[!, var]) <:Union{Missing, Number})
            non_num_str = non_num_str*string(var)*", "
        end
    end
    if non_num_str != ""
        throw("The following columns should be numeric but are not: "*non_num_str[1:end-2])
    end

end

function error_throw_fnc_sec_stage(data::DataFrame, 
    flexible_input::Symbol, 
    fixed_inputs::Union{Symbol,Array{Symbol}}, 
    id::Symbol, 
    time::Symbol,
    opts::Dict)

    # Check if variables are actually in DataFrame
    missing_str = string()
    non_num_str = string()

    all_var_symbols = vec(hcat(fixed_inputs..., flexible_input, id, time))

    for var in all_var_symbols
        if string(var) ∉ names(data)
            missing_str = missing_str*string(var)*", "
        end
    end
    if missing_str != ""
        throw("The following columns are not in the specified dataframe: "*missing_str[1:end-2])
    end

    # Check if inputs and output are numeric
    for var in setdiff(all_var_symbols, vec(hcat(id, time)))
        if !(eltype(data[!, var]) <:Union{Missing, Number})
            non_num_str = non_num_str*string(var)*", "
        end
    end
    if non_num_str != ""
        throw("The following columns should be numeric but are not: "*non_num_str[1:end-2])
    end

end

## Auxiliary function to fill up options that were not given in opts dictionary
function opts_filler(opts::Dict)

    # Define new opts dictionary because opts can have a too narrow type if the user did only specify a specific subset of options
    new_opts::Dict{String, Any} = opts
    
    if "fes_print_starting_values" ∉ keys(new_opts)
        # opts["fes_print_starting_values"] = false
        new_opts["fes_print_starting_values"] = false
    end
    if "fes_print_results" ∉ keys(new_opts)
        new_opts["fes_print_results"] = false
    end
    if "fes_method" ∉ keys(new_opts)
        new_opts["fes_method"] = "NLLS"
    end
    if "fes_optimizer" ∉ keys(new_opts)
        new_opts["fes_optimizer"] = NelderMead()
    end
    if "fes_optimizer_options" ∉ keys(new_opts)
        new_opts["fes_optimizer_options"] = Optim.Options(iterations = 20000,
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
    if "ses_print_starting_values" ∉ keys(new_opts)
        new_opts["ses_print_starting_values"] = false
    end
    if "ses_optimizer" ∉ keys(new_opts)
        new_opts["ses_optimizer"] = NelderMead()
    end
    if "ses_optimizer_options" ∉ keys(new_opts)
        new_opts["ses_optimizer_options"] = Optim.Options(iterations = 20000,
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
    if "ses_print_starting_values" ∉ keys(new_opts) 
        new_opts["ses_print_starting_values"] = false
    end
    if "ses_print_results" ∉ keys(new_opts)
        new_opts["ses_print_results"] = false
    end

    if "print_results" ∉ keys(new_opts)
        opts["print_results"] = true
    end

    if "maxboottries" ∉ keys(new_opts)
        new_opts["maxboottries"] = 10
    end

    if "called_from_bootstrapping" ∉ keys(new_opts)
        opts["called_from_bootstrapping"] = false
    end
    
    return new_opts
end

## Function that checks a string to only contain defined substrings
function check_str_only_def_substr(s::String, strings_to_check::Union{Vector,String,Char})
    parts = split(s, '⋅')

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

## Function that iterates over an array of strings and checks if it only contains defined substrings (or "⋅")
function check_array_string_only_substrings(;s_vec::Vector{String}, strings_to_check::Union{Vector,String,Char})
    check_res = falses(length(s_vec))
    j = 1
    for stri in s_vec
        res = check_str_only_def_substr(stri, [strings_to_check..., '⋅'])
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
        for sym in all_var_symbols # Iterate over all symbols of the polynomials
            parts = split(string(sym), '⋅')
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
function startvalues(;data::DataFrame, Y_var::Symbol, X_vars::Union{Symbol,Array{Symbol}}, user_start_vals::Vector{<:Union{Missing, Number}}, stage::String, opts::Dict)
    # Check if user-provided starting values are in the right form
    if !(all(ismissing.(user_start_vals)) && length(user_start_vals) == 1)
        # If user-provided not the right number of starting values
        if length(user_start_vals) != length(X_vars)
            throw(uppercasefirst(stage)*" starting values: You did not provide the correct number of starting values. They need to match the number of terms in the polynomial series. You can also leave them unspecified to let the program choose starting values.")
        end
    end

    if any(ismissing.(user_start_vals)) || size(user_start_vals) == (1,) # If user did not specify starting values

        # Define matrices for OLS
        Y = data[:, Y_var]
        X = hcat(data.constant, Matrix(data[:,X_vars]))

        # Calculate OLS
        startvals = vec(inv(X'*X)*X'*Y)

        if stage == "first stage" # GNR make a correction to the constant in their first stage (this is probably because the constant is so strongly negative for their example that the ln(Xγ) returns NaNs because of negative Xγ-values)
            startvals[1] = 0.1 - minimum((X*startvals) .- startvals[1])
        end
    else # If user specified starting starting_values

        # Check dimensions and throw an error if they do not match
        if size(X_vars)[1] + 1 != size(user_start_vals)[1] # Add constant
            println(size(X_vars))
            println(size(user_start_vals))
            throw("You specified either too many or not enough "*stage*" starting values!")

        else
            startvals = user_start_vals
        end
    end

    # Print results if specified
    if ((opts["fes_print_starting_values"] == true && stage == "first stage") || (opts["ses_print_starting_values"] == true && stage == "second stage")) && opts["called_from_bootstrapping"] == false
        
        print_tab = hcat(vcat([:constant], X_vars), startvals)
        if stage == "second stage"
            print_tab = hcat(vcat([:constant], X_vars), startvals)[begin+1:end, :] # Because we don't need the constant in the second stage
        end
        header = (["Variable", "Value"])

        # Make first character uppercase if it is lowercase
        if islowercase(stage[1]) == true
            stage = uppercasefirst(stage)
        end

        println(stage*" starting values:")
        pretty_table(print_tab, header = header, formatters =  ft_printf("%5.5f"), limit_printing = false)
    end

    # Return result
    return startvals
end

## Polynomial series generating function
function polynom_series!(;data::DataFrame, var_names::Union{Vector{Symbol},Symbol}, degree::Int)

    # Check if user put in an invalid degree
    if degree < 1
        throw(error("Polynomial series degree must be at least 1"))
    end

    # Convert to a vector if only one symbol put in to make my life easier
    if typeof(var_names) == Symbol
        var_names = [var_names]
    end

    # Preallocate ooutput
    poly_var_names = Set(Symbol[])

    # Recursively define polynomail calc fnc
    function generate_polynomials(data, var_names, degree, prefix, idx)

        # Stop if degree is 0
        if degree == 0
            return
        end

        # Calculate poynomials
        for i in idx:length(var_names)
            new_prefix = prefix == "" ? string(var_names[i]) : string(prefix, '⋅', var_names[i])
            varn = Symbol(new_prefix)
            push!(poly_var_names, varn)
            if varn ∈ names(data)
                data = select(data, Not(varn))
            end

            if prefix == ""
                data[:, varn] = data[:, var_names[i]]
            else
                prev_var = Symbol(prefix)
                data[:, varn] = data[:, prev_var] .* data[:, var_names[i]]
            end

            generate_polynomials(data, var_names, degree - 1, new_prefix, i)
        end
    end

    # Get vector of o
    for var in var_names
        push!(poly_var_names, var)
    end

    # Run the just-defined fnc
    generate_polynomials(data, var_names, degree, "", 1)
    
    return collect(poly_var_names)
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
    if any(in(Symbol.(names(data))), Symbol.(lag_prefix, variable))
        if force == true
            data = data[!, Not(filter(in(Symbol.(names(data))), Symbol.(lag_prefix, variable)))]
        else
            throw("Specified name for lag of variable already present in specified dataframe. Either set force = true, choose difference lag variable name, or rename the column.")
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

    if any(in(Symbol.(names(data))), Symbol.(lag_prefix, variable))
        if force == true
            select!(data, Not(filter(in(Symbol.(names(data))), Symbol.(lag_prefix, variable))))
            # df2 = df2[!, Not(Symbol.(lag_prefix, variable))]
        else
            throw("Specified name for lag of variable already present in specified dataframe. Either set force = true, choose difference lag variable name, or rename the column.")
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

## Function that returns the superscript of the corresponding input character (b/c Julia does not have a simple function for that)
# Define a dictionary for superscript characters
const superscript_map = Dict(
    '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
    '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
    'a' => 'ᵃ', 'b' => 'ᵇ', 'c' => 'ᶜ', 'd' => 'ᵈ', 'e' => 'ᵉ',
    'f' => 'ᶠ', 'g' => 'ᵍ', 'h' => 'ʰ', 'i' => 'ⁱ', 'j' => 'ʲ',
    'k' => 'ᵏ', 'l' => 'ˡ', 'm' => 'ᵐ', 'n' => 'ⁿ', 'o' => 'ᵒ',
    'p' => 'ᵖ', 'r' => 'ʳ', 's' => 'ˢ', 't' => 'ᵗ', 'u' => 'ᵘ',
    'v' => 'ᵛ', 'w' => 'ʷ', 'x' => 'ˣ', 'y' => 'ʸ', 'z' => 'ᶻ',
    '+' => '⁺', '-' => '⁻', '=' => '⁼', '(' => '⁽', ')' => '⁾'
)

function superscript_this!(c::String) # Need to use a string as input because I don't understand Chars in Julia. Char(5) returns a different unicode than string(5). And the superscript of Char(5) does not  work
    # Return the superscript character if it exists in the map, else return the original character
    return get(superscript_map, c[1], c[1])
end