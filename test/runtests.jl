# Load required dependencies (use 'using Pkg' and Pkg.add("DataFrames"), Pkg.add("GNRProdEst") or Pkg.add("CSV") if you do not have any of these packages in your environment.)
using GNRProdEst, DataFrames, CSV, Test, Optim


@testset "GNRProdEst.jl" begin

    # Read in replication data. @__DIR__ is this file's directory, so the data is found no
    # matter what the working directory is. Guessing from pwd() only worked under Pkg.test
    # (which cds into test/) and left rep_data undefined when run from the package root.
    rep_data = CSV.read(joinpath(@__DIR__, "GNR_data_500.csv"), DataFrame)

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
                                                            boot_reps = 10,
                                                            opts = opts
                                                    );

    @test round.(gnr_fes_res["γ"], digits = 5) == [0.65239, -0.00112, 0.00476, -0.00146, 0.00111, -5.0e-5, -0.0005, 0.00123, -0.00033, -0.02119]

    @test round.(gnr_fes_res["γ_dash"], digits = 5) == [0.6759, -0.00116, 0.00494, -0.00152, 0.00115, -5.0e-5, -0.00052, 0.00128, -0.00035, -0.02195]

    @test round.(gnr_fes_res["γ_flex"], digits = 5) == [0.65239, -0.00037, 0.00159, -0.00146, 0.00055, -1.0e-5, -0.0005, 0.00062, -0.00033, -0.01059]

    @test gnr_fes_res["E"] ≈ 1.03604299 rtol = 1e-8

    @test Set(gnr_fes_res["polynom_series"]) == Set(Symbol.(["k"; "k⋅k"; "k⋅k⋅k"; "k⋅k⋅i"; "k⋅i"; "k⋅i⋅i"; "i"; "i⋅i"; "i⋅i⋅i"]))

    @test Set(gnr_fes_res["all_inputs"]) == Set([:k; :i])

    @test gnr_fes_res["fixed_inputs"] == [:k]

    @test gnr_fes_res["flexible_input"] == :i

    @test gnr_fes_res["share_degree"] == 3

    @test size(gnr_fes_res["fes_optim_estimates"]) == (10, 4)

    @test typeof(gnr_fes_res["fes_optim_estimates"]) == DataFrame

    @test names(gnr_fes_res["fes_optim_estimates"]) == ["Variable"; "γ"; "γ_dash"; "γ_flex"]

    @test Optim.converged(gnr_fes_res["fes_optim_results"])

    @test round.(gnr_ses_res["α"], digits = 5) == [0.38825, -0.02485, 0.00247]

    @test round.(vec(gnr_ses_res["δ"]), digits = 4) == [0.168, 0.7691, 0.0654, -0.0398] # 4 digits: δ[4] differs in the 5th digit between a plain run and one under --check-bounds=yes (which Pkg.test uses)
    
    @test Optim.converged(gnr_ses_res["gmm_optim_results"]) 

    @test Set(gnr_ses_res["polynom_fixed"]) == Set(Symbol.(["k"; "k⋅k"; "k⋅k⋅k"]))

    @test gnr_ses_res["fixed_inputs"] == [:k]

    @test gnr_ses_res["flexible_input"] == :i

    @test gnr_ses_res["int_const_series_degree"] == 3

    @test gnr_ses_res["lm_tfp_degree"] == 3

    @test Set(gnr_ses_res["all_inputs"]) == Set([:k; :i])
end
