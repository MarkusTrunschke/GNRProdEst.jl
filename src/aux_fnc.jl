## Function that converts inputs into correct types for my program
"""
    GNR_input_cleaner!(; fixed_inputs, flexible_input, fes_starting_values=[missing], ses_starting_values=[missing]) -> Tuple

Bring user-supplied arguments into the shapes the internals expect. Wraps a single
`Symbol` into a one-element vector, flattens matrices of symbols into vectors, and
converts integer starting values to floats because the optimizers reject integers.

# Keyword Arguments
- `fixed_inputs::Union{Array{Symbol},Symbol}`: Fixed input variable(s)
- `flexible_input::Union{Array{Symbol},Symbol}`: Flexible input variable
- `fes_starting_values::Union{Vector{<:Number},Vector{Missing}}=[missing]`: First stage starting values
- `ses_starting_values::Union{Vector{<:Number},Vector{Missing}}=[missing]`: Second stage starting values

# Returns
- `Tuple`: `(fixed_inputs, flexible_input, fes_starting_values, ses_starting_values)`, cleaned
"""
function GNR_input_cleaner!(;fixed_inputs::Union{Array{Symbol},Symbol}, flexible_input::Union{Array{Symbol},Symbol}, fes_starting_values::Union{Vector{<:Number},Vector{Missing}} = vec([missing]), ses_starting_values::Union{Vector{<:Number},Vector{Missing}} = vec([missing]))
    # Convert fixed input into vector of symbols if an array was given
    if fixed_inputs isa AbstractArray # isa rather than typeof(...) == Array: typeof returns a concrete type, so the old test was never true and only Matrix{Symbol} was flattened. It also keeps vec, which has no method for a Symbol, provably reachable only for arrays
        fixed_inputs = vec(fixed_inputs)
    end

    # Convert fixed input to a vector if user put in a symbol. (Makes working with it in the program easier) 
    if fixed_inputs isa Symbol
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
"""
    error_throw_fnc(data, output, flexible_input, fixed_inputs, ln_share_flex_y_var, id, time, opts)

Validate the arguments of the combined estimation routine. Throws if a requested column
is missing from `data` or if a column that has to be numeric is not. Returns `nothing`.

# Arguments
- `data::DataFrame`: Dataset to check
- `output::Symbol`: Output variable column name
- `flexible_input::Symbol`: Flexible input variable
- `fixed_inputs::Union{Symbol,Array{Symbol}}`: Fixed input variable(s)
- `ln_share_flex_y_var::Symbol`: Log flexible input share of output
- `id::Symbol`: Firm identifier
- `time::Symbol`: Time identifier
- `opts::Dict`: Further options
"""
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

"""
    error_throw_fnc_first_stage(data, output, flexible_input, fixed_inputs, ln_share_flex_y_var, opts)

Validate the arguments of the first stage. Same checks as `error_throw_fnc`, without the
panel identifiers, which the first stage does not need. The share variable is only checked
when the user supplied one. Returns `nothing`.

# Arguments
- `data::DataFrame`: Dataset to check
- `output::Symbol`: Output variable column name
- `flexible_input::Symbol`: Flexible input variable
- `fixed_inputs::Union{Symbol,Array{Symbol}}`: Fixed input variable(s)
- `ln_share_flex_y_var::Symbol`: Log flexible input share of output, or `:NotDefinedByUser`
- `opts::Dict`: Further options
"""
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

"""
    error_throw_fnc_sec_stage(data, flexible_input, fixed_inputs, id, time, opts)

Validate the arguments of the second stage. Throws if a requested column is missing from
`data` or is not numeric. Returns `nothing`.

# Arguments
- `data::DataFrame`: Dataset to check
- `flexible_input::Symbol`: Flexible input variable
- `fixed_inputs::Union{Symbol,Array{Symbol}}`: Fixed input variable(s)
- `id::Symbol`: Firm identifier
- `time::Symbol`: Time identifier
- `opts::Dict`: Further options
"""
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
"""
    opts_filler(opts::Dict) -> Dict{String,Any}

Fill in every option the package needs that the user did not specify. Returns a new
`Dict{String,Any}`, because a dictionary built from a few options only is typically typed
too narrowly to hold the optimizer objects and strings added here.

Defaults set are the print flags (all `false` except `print_results`), `fes_method`
(`"NLLS"`), the first and second stage optimizers (`NelderMead()`) and their
`Optim.Options`, `maxboottries` (`10`), and `called_from_bootstrapping` (`false`).

# Arguments
- `opts::Dict`: Options given by the user

# Returns
- `Dict{String,Any}`: The user's options with all missing entries filled in
"""
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
                                                      f_reltol = 1e-9, # Was f_tol, deprecated by Optim in favour of f_reltol / f_abstol. Optim mapped f_tol onto f_reltol with the same value, so behaviour is unchanged
                                                      x_abstol = 1e-12, # Was x_tol, deprecated by Optim in favour of x_abstol / x_reltol. Optim mapped x_tol onto x_abstol with the same value, so behaviour is unchanged
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
                                                      f_reltol = 1e-9, # Was f_tol, deprecated by Optim in favour of f_reltol / f_abstol. Optim mapped f_tol onto f_reltol with the same value, so behaviour is unchanged
                                                      x_abstol = 1e-12, # Was x_tol, deprecated by Optim in favour of x_abstol / x_reltol. Optim mapped x_tol onto x_abstol with the same value, so behaviour is unchanged
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
        new_opts["print_results"] = true # Must be new_opts, not opts, for the same reason as called_from_bootstrapping below
    end

    if "maxboottries" ∉ keys(new_opts)
        new_opts["maxboottries"] = 10
    end

    if "called_from_bootstrapping" ∉ keys(new_opts)
        new_opts["called_from_bootstrapping"] = false # Must be new_opts, not opts: new_opts is a copy whenever the caller passes a narrower Dict (e.g. Dict{String,Bool}), and new_opts is what gets returned
    end
    
    return new_opts
end

## Function that checks a string to only contain defined substrings
"""
    check_str_only_def_substr(s::String, strings_to_check) -> Bool

Test whether a polynomial name is built exclusively from the given variables. Splits `s`
on the internal separator `⋅` and returns `true` only if every part matches one of
`strings_to_check`.

# Arguments
- `s::String`: Polynomial name, e.g. `"k⋅k⋅i"`
- `strings_to_check::Union{Vector,String,Char}`: Variable name(s) that are allowed to appear

# Returns
- `Bool`: `true` if all parts of `s` are in `strings_to_check`
"""
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
"""
    check_array_string_only_substrings(; s_vec, strings_to_check) -> Matrix

Apply `check_str_only_def_substr` to a vector of polynomial names. Used to pick the pure
fixed input polynomials out of the full first stage series.

# Keyword Arguments
- `s_vec::Vector{String}`: Polynomial names to check
- `strings_to_check::Union{Vector,String,Char}`: Variable name(s) that are allowed to appear

# Returns
- `Matrix`: Two columns, the names and a `Bool` per name
"""
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
"""
    get_input_degree(input, all_var_symbols) -> Matrix

Determine how often each input appears in each polynomial of a series, i.e. the degree of
that input in that term. Counts occurrences in the internal naming scheme, in which the
components of a polynomial are separated by `⋅`.

# Arguments
- `input::Union{Symbol,Vector{Symbol}}`: Input(s) whose degree is counted
- `all_var_symbols::Union{Symbol,Vector{Symbol}}`: Polynomial series names

# Returns
- `Matrix`: First row holds the polynomial names, row `1 + i` the degree of `input[i]` in each
"""
function get_input_degree(input::Union{Symbol,Vector{Symbol}}, all_var_symbols::Union{Symbol,Vector{Symbol}})
    
    # Clean inputs. The annotated locals keep the types provable for the length calls below,
    # neither of which has a method for a lone Symbol
    input_vec::Vector{Symbol} = input isa Symbol ? [input] : input
    all_vars::Vector{Symbol} = all_var_symbols isa Symbol ? [all_var_symbols] : all_var_symbols

    # Initalize array
    input_degree_mat = Array{Union{Symbol,<:Number}}(undef, length(input_vec) + 1, length(all_vars))# Array{Union{Symbol,Int}}[]
    
    # Give first row all variable symbols
    input_degree_mat[1,:] = all_vars

    i = 2
    for inp in input_vec # Loop over all inputs
        j = 1
        for sym in all_vars # Iterate over all symbols of the polynomials
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
"""
    startvalues(; data, Y_var, X_vars, user_start_vals, stage, opts) -> Vector

Return starting values for one of the two estimation stages. If the user supplied values
they are passed through unchanged, otherwise they are calculated from an OLS regression of
`Y_var` on `X_vars`, following the GNR replication code. In the first stage the constant is
corrected upwards so that the fitted values entering the logarithm are positive.

# Keyword Arguments
- `data::DataFrame`: Dataset
- `Y_var::Symbol`: Dependent variable of the auxiliary regression
- `X_vars::Vector{Symbol}`: Independent variables of the auxiliary regression
- `user_start_vals::Vector{<:Union{Missing, Number}}`: User-supplied starting values, or `[missing]`
- `stage::String`: `"first stage"` or `"second stage"`
- `opts::Dict`: Further options, used to decide whether to print the values

# Returns
- `Vector`: Starting values, including the constant
"""
function startvalues(;data::DataFrame, Y_var::Symbol, X_vars::Vector{Symbol}, user_start_vals::Vector{<:Union{Missing, Number}}, stage::String, opts::Dict)
    # Check if user-provided starting values are in the right form
    if !(all(ismissing.(user_start_vals)) && length(user_start_vals) == 1)
        # If user-provided not the right number of starting values
        if length(user_start_vals) != length(X_vars) + 1 # + 1 for the constant. Without it this check and the one further down contradict each other, so no number of starting values could ever pass both
            throw(uppercasefirst(stage)*" starting values: You did not provide the correct number of starting values. They need to match the number of terms in the polynomial series plus one for the constant. You can also leave them unspecified to let the program choose starting values.")
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
        if length(X_vars) + 1 != length(user_start_vals) # Add constant
            println(length(X_vars))
            println(length(user_start_vals))
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
        pretty_table(print_tab, column_labels = header, formatters = [fmt__printf("%5.5f")], limit_printing = false)
    end

    # Return result
    return startvals
end

## Polynomial series generating function
"""
    polynom_series!(; data, var_names, degree) -> Vector{Symbol}

Add all polynomials and interactions of `var_names` up to `degree` to `data` in-place. The
new columns are named by joining the components with `⋅`, so that the degree of each input
in each term can be recovered later by `get_input_degree`.

# Keyword Arguments
- `data::DataFrame`: Data frame to mutate
- `var_names::Union{Vector{Symbol},Symbol}`: Variable(s) to expand
- `degree::Int`: Highest polynomial degree, must be at least 1

# Returns
- `Vector{Symbol}`: Names of the generated columns
"""
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
"""
    polynomial_fnc_fast!(poly_mat, degree; par_cal=false) -> poly_mat

In-place computation of polynomial terms. Fills columns 2 through `degree` of `poly_mat` with powers 2 through `degree` of column 1. Preallocated matrix avoids allocations during repeated calls.

# Arguments
- `poly_mat::AbstractArray{<:Number}`: Matrix with base values in column 1
- `degree::Int`: Highest polynomial degree to compute

# Keyword Arguments
- `par_cal::Bool=false`: Use `Threads.@threads` for parallel computation
"""
function polynomial_fnc_fast!(poly_mat::AbstractArray{<:Number}, degree::Int; par_cal::Bool = false) # AbstractArray covers both Array and the SubArray produced by @view at the call sites, and unlike Union{Array{<:Number},SubArray{<:Number}} it is resolvable by static analysers
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
"""
    panel_lag(; data, id, time, variable, lag_prefix="lag_", lags=1, drop_missings=false, force=false) -> DataFrame

Non-mutating version of `panel_lag!`. Copies `data` before lagging, so the caller's data
frame is left untouched, and returns the copy.

# Keyword Arguments
- `data::DataFrame`: Data frame to lag
- `id::Symbol`: Panel identifier column
- `time::Symbol`: Time variable for ordering
- `variable::Union{Symbol,Vector{Symbol}}`: Column(s) to lag
- `lag_prefix::String="lag_"`: Prefix for lag column names
- `lags::Int=1`: Lag distance
- `drop_missings::Bool=false`: Drop rows with missing lags
- `force::Bool=false`: Remove existing lag columns if present

# Returns
- `DataFrame`: Copy of `data` with the lag columns added
"""
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
"""
    panel_lag!(; data, id, time, variable, lag_prefix="lag_", lags=1, drop_missings=false, force=false) -> DataFrame

Compute panel lags in-place using `ShiftedArrays.lag` within groups. Sorts by `id` and `time`, creates lag columns with specified prefix, and validates time gaps.

# Keyword Arguments
- `data::DataFrame`: Data frame to mutate
- `id::Symbol`: Panel identifier column
- `time::Symbol`: Time variable for ordering
- `variable::Union{Array{Symbol},Symbol}`: Column(s) to lag
- `lag_prefix::String="lag_"`: Prefix for lag column names
- `lags::Int=1`: Lag distance
- `drop_missings::Bool=false`: Drop rows with missing lags
- `force::Bool=false`: Remove existing lag columns if present
"""
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

"""
    lagging_that_panel!(; data, id, time, variable, lag_prefix="lag_", lags=1, drop_missings=false) -> DataFrame

Internal helper for `panel_lag!`. Computes lags using `ShiftedArrays.lag` within groups, validates time gaps (sets to `missing` if gap ≠ `lags`), joins back to data, and renames columns.
"""
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

"""
    superscript_this!(c::String) -> Char

Convert first character of string to its Unicode superscript equivalent using `superscript_map`. Returns original character if no superscript exists.
"""
function superscript_this!(c::String) # Need to use a string as input because I don't understand Chars in Julia. Char(5) returns a different unicode than string(5). And the superscript of Char(5) does not  work
    # Return the superscript character if it exists in the map, else return the original character
    return get(superscript_map, c[1], c[1])
end