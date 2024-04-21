## General estimation command
function GNRProd(data::DataFrame, output::Symbol, flex_inputs::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}})
    
    ## Run data preparation to ensure inputs are working with internals
    est_df = prep_data(data; output = output, flex_inputs = flex_inputs, fixed_inputs = fixed_inputs)

    ## Run first stage estimation
    sol = GNRFirstStage(est_df; output = output, flex_inputs = flex_inputs, fixed_inputs = fixed_inputs, share = :share_m_y)

    return sol
end
## First stage estimation
function GNRFirstStage(data; output::Symbol, flex_inputs::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, share::Symbol )

    # Get starting values for NLLS (Follow replication code of GNR for now)
    γ_start = FS_startvalues(data)


    # Test function
    function f(u, p)
        u .* u .- p
    end
    u0 = [1.0, 1.0]
    p = 2.0
    prob = NonlinearProblem(f, u0, p)
    sol = solve(prob, NewtonRaphson())
    print(sol)

    return sol
end

function FS_startvalues(data, share_var, flex_inputs, fixed_inputs)
    
end
