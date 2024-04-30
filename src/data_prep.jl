function prep_data(data::DataFrame; output::Symbol, flex_input::Symbol, fixed_inputs::Union{Symbol,Array{Symbol}})
    
    ## Select necessary variables from data frame
    # Flatten all arguments. Some might be arrays of symbols and others might just be symbols. The following iterates over all sublists and flattens them
    all_var_symbols = [x for sublist in [output, flex_input, fixed_inputs] for x in (sublist isa Vector ? sublist : [sublist])]

    all_input_symbols = [x for sublist in [flex_input, fixed_inputs] for x in (sublist isa Vector ? sublist : [sublist])]
    
    # Select the data
    est_df = data[:,all_var_symbols]

    # Add a constant to the data
    est_df.constant = ones(size(data)[1])

    # Drop missings
    est_df = dropmissing(est_df) # Drop missings

    # Calculate share of intermediate input on revenue variable
    est_df[!, :share_flex_y] = select(est_df, flex_input)[:,1] ./ select(est_df, output)[:,1] # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...
    est_df[!, :ln_share_flex_y] = log.(select(est_df, flex_input)[:,1] ./ select(est_df, output)[:,1]) # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...

    # Return the data
    return est_df, all_var_symbols, all_input_symbols
end