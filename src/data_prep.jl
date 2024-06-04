function prep_data!(data::DataFrame; output::Symbol, flexible_input::Symbol, fixed_inputs::Union{Symbol,Array{Symbol}}, ln_share_flex_y_var::Symbol, id::Symbol, time::Symbol)
    
    ## Select necessary variables from data frame
    # Flatten all arguments. Some might be arrays of symbols and others might just be symbols. The following iterates over all sublists and flattens them
    all_var_symbols = [x for sublist in [output, flexible_input, fixed_inputs, ln_share_flex_y_var, id, time] for x in (sublist isa Vector ? sublist : [sublist])]
    all_input_symbols = [x for sublist in [flexible_input, fixed_inputs] for x in (sublist isa Vector ? sublist : [sublist])]
    
    # Select the data
    select!(data, all_var_symbols)

    # Add a constant to the data
    data.constant = ones(size(data)[1])

    # Drop missings
    dropmissing!(data)

    # Calculate share of intermediate input on revenue variable
    if ln_share_flex_y_var == :NotDefinedByUser
        data[!, :ln_share_flex_y] = select(data, flexible_input)[:,1] .- select(data, output)[:,1] # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...
        ln_share_flex_y_var = :ln_share_flex_y
    end

    # Return the data
    return data, all_input_symbols, ln_share_flex_y_var
end