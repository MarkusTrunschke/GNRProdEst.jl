function prep_data(data::DataFrame; output::Symbol, flex_inputs::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, intermediate_input::Symbol)
    
    ## Select necessary variables from data frame
    # Flatten all arguments. Some might be arrays of symbols and others might just be symbols. The following iterates over all sublists and flattens them
    all_var_symbols = [x for sublist in [output, flex_inputs, fixed_inputs] for x in (sublist isa Vector ? sublist : [sublist])]
    
    # Select the data
    est_df = data[:,all_var_symbols]

    # Drop missings
    est_df = dropmissing(est_df) # Drop missings

    # Calculate share of intermediate input on revenue variable
    est_df[!, :ln_share_m_y] = log.(select(est_df, intermediate_input)[:,1] ./ select(est_df, output)[:,1]) # Need to do the [:,1] b/c need to convert it to a vector before adding in to a column...

    # Return the data
    return est_df
end