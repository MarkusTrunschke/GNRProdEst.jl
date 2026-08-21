"""
    prep_data!(data; output, flexible_input, fixed_inputs, ln_share_flex_y, id, time) -> Tuple

Prepare a dataset for estimation in-place. Reduces `data` to the columns the estimation
needs, adds a constant column, drops rows with missing values, and derives the log
flexible input share of output when the user did not supply it, as the difference between
the log flexible input and log output.

# Arguments
- `data::DataFrame`: Data frame to mutate

# Keyword Arguments
- `output::Symbol`: Log output variable
- `flexible_input::Symbol`: Log flexible input variable
- `fixed_inputs::Union{Symbol,Array{Symbol}}`: Log fixed input variable(s)
- `ln_share_flex_y::Symbol`: Log flexible input share of output, or `:NotDefinedByUser` to derive it
- `id::Symbol`: Firm identifier
- `time::Symbol`: Time identifier

# Returns
- `Tuple`: `(data, all_input_symbols, ln_share_flex_y)`, where the last element is the name of
  the share column actually used
"""
function prep_data!(data::DataFrame; output::Symbol, flexible_input::Symbol, fixed_inputs::Union{Symbol,Array{Symbol}}, ln_share_flex_y::Symbol, id::Symbol, time::Symbol)
    
    ## Select necessary variables from data frame
    # Flatten all arguments. Some might be arrays of symbols and others might just be symbols. The following iterates over all sublists and flattens them
    # all_var_symbols = [x for sublist in [output, flexible_input, fixed_inputs, ln_share_flex_y_var, id, time] for x in (sublist isa Vector ? sublist : [sublist])]
    # all_input_symbols = [x for sublist in [flexible_input, fixed_inputs] for x in (sublist isa Vector ? sublist : [sublist])]
    all_var_symbols = Array{Symbol}(undef,0)
    if ln_share_flex_y != :NotDefinedByUser
        all_var_symbols = vec(hcat(output, flexible_input, fixed_inputs... , ln_share_flex_y, id, time))
    else
        all_var_symbols = vec(hcat(output, flexible_input, fixed_inputs... , id, time))
    end

    all_input_symbols = vec(hcat(flexible_input, fixed_inputs...))
    # Select the data
    select!(data, all_var_symbols)

    # Add a constant to the data
    data.constant = ones(size(data)[1])

    # Drop missings
    dropmissing!(data)

    # Calculate share of intermediate input on revenue variable
    if ln_share_flex_y == :NotDefinedByUser
        data[!, :ln_share_flex_y] = select(data, flexible_input)[:,1] .- select(data, output)[:,1] # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...
        ln_share_flex_y = :ln_share_flex_y
    end

    # Return the data
    return data, all_input_symbols, ln_share_flex_y
end