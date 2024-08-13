# Load required dependencies (use 'using Pkg' and Pkg.add("DataFrames"), Pkg.add("GNRProdEst") or Pkg.add("CSV") if you do not have any of these packages in your environment.)
using GNRProdEst, DataFrames, CSV, Test, Optim


@testset "GNRProdEst.jl" begin

    # Read in replication data
    if last(pwd(), 4) == "test"
        rep_data = CSV.read(joinpath(pwd(),"GNR_data_500.csv"), DataFrame)
    elseif last(pwd(), 12) == "NRProdEst.jl"
        rep_data = CSV.read(joinpath(pwd(),"test","GNR_data_500.csv"), DataFrame)
    end

    # Define some options to print results
    opts = Dict("ses_optimizer_options" => Optim.Options(f_tol = 1e-12,
                                                         x_tol = 1e-12,
                                                         g_tol = 1e-12),
                "fes_print_results" => true,
                "fes_print_starting_values" => true
                )

    # Run both estimation stages at the same time
    gnr_fes_res, gnr_ses_res = GNRProdEst.gnrprodest(data = rep_data, 
                                                            output = :yg, 
                                                            flexible_input = :i, 
                                                            fixed_inputs = :k, 
                                                            ln_share_flex_y = :si, 
                                                            id = :id, 
                                                            time = :time,
                                                            opts = opts
                                                    );

    @test gnr_fes_res["γ"] ≈  [0.6523861637358107
                              -0.0011171888341931346
                               0.004764634033341126
                              -0.0014626656866083623
                               0.001107644701812702
                              -4.5080251351363865e-5
                              -0.0004990739219499631
                               0.0012333049857546326
                              -0.000334609441120487
                              -0.021189821854455733]

    @test gnr_fes_res["γ_dash"] ≈ [ 0.6759001109073727
                                   -0.001157455658795065
                                    0.0049363656842857945
                                   -0.0015153845295213049
                                    0.0011475675273584498
                                   -4.6705078344454365e-5
                                   -0.0005170620377129271
                                    0.0012777569835030124
                                   -0.0003466697654482715
                                   -0.021953566365539898]

    @test gnr_fes_res["γ_flex"] ≈ [ 0.6523861637358107
                                   -0.0003723962780643782
                                    0.001588211344447042
                                   -0.0014626656866083623
                                    0.000553822350906351
                                   -1.1270062837840966e-5
                                   -0.0004990739219499631
                                    0.0006166524928773163
                                   -0.000334609441120487
                                   -0.010594910927227866]

    @test gnr_fes_res["E"] == 1.0360429887674383

    @test Set(gnr_fes_res["polynom_series"]) == Set(Symbol.(["k"; "k⋅k"; "k⋅k⋅k"; "k⋅k⋅i"; "k⋅i"; "k⋅i⋅i"; "i"; "i⋅i"; "i⋅i⋅i"]))

    @test Set(gnr_fes_res["all_inputs"]) == Set([:k; :i])

    @test gnr_fes_res["fixed_inputs"] == [:k]

    @test gnr_fes_res["flexible_input"] == :i

    @test gnr_fes_res["share_degree"] == 3

    @test size(gnr_fes_res["fes_optim_estimates"]) == (10, 4)

    @test typeof(gnr_fes_res["fes_optim_estimates"]) == DataFrame

    @test names(gnr_fes_res["fes_optim_estimates"]) == ["Variable"; "γ"; "γ_dash"; "γ_flex"]

    @test Optim.converged(gnr_fes_res["fes_optim_results"])

    @test gnr_ses_res["α"] == [0.3871439515587196, -0.02456417904454124, 0.0024446787683550806]

    @test gnr_ses_res["δ"] == [0.16830835612691358; 0.7688330072496795; 0.06571422537798424; -0.03985539721866669;;]
    
    @test Optim.converged(gnr_ses_res["gmm_optim_results"]) 

    @test Set(gnr_ses_res["polynom_fixed"]) == Set(Symbol.(["k"; "k⋅k"; "k⋅k⋅k"]))

    @test gnr_ses_res["fixed_inputs"] == [:k]

    @test gnr_ses_res["flexible_input"] == :i

    @test gnr_ses_res["int_const_series_degree"] == 3

    @test gnr_ses_res["lm_tfp_degree"] == 3

    @test Set(gnr_ses_res["all_inputs"]) == Set([:k; :i])
end