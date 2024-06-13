@testset "GNRProdEst.jl" begin

    # Load required dependencies (use 'using Pkg' and Pkg.add("DataFrames"), Pkg.add("GNRProdEst") or Pkg.add("CSV") if you do not have any of these packages in your environment.)
    using GNRProdEst, DataFrames, CSV, Test, Optim

    # Read in replication data
    rep_data = CSV.read("test/GNR_data_500.csv", DataFrame)

    # Define some options to print results
    opts = Dict("fes_print_results" => true, "ses_print_results" => true)

    # Run both estimation stages at the same time
    gnr_fes_res, gnr_ses_res = GNRProdEst.gnrprodest(data = rep_data, 
                                    output = :yg, 
                                    flexible_input = :i, 
                                    fixed_inputs = :k, 
                                    ln_share_flex_y = :si, 
                                    id = :id, 
                                    time = :time,
                                    opts = Dict("fes_print_results" => false,
                                                "ses_print_results" => false));

    @test gnr_fes_res["γ"] == [0.6523861561619183
                               -4.50793426512781e-5
                               -0.0004990736355941706
                               -0.0014626620369507918
                                0.001233304692081355
                               -0.00033461019716388846
                                0.0011076470199486048
                                0.004764633570968514
                               -0.021189821397540238
                               -0.0011171913074503068]

    @test gnr_fes_res["γ_dash"] == [ 0.675900102838662
                                    -4.6704136876773166e-5
                                    -0.0005170617408663145
                                    -0.0015153807478218147
                                     0.0012777566788255083
                                    -0.00034667054862795837
                                     0.0011475699286702626
                                     0.00493636520362776
                                    -0.021953565884950572
                                    -0.0011574582208159364]

    @test gnr_fes_res["γ_flex"] == [0.6523861561619183
                                    -1.1269835662819525e-5
                                    -0.0004990736355941706
                                    -0.0014626620369507918
                                     0.0006166523460406775
                                    -0.00033461019716388846
                                     0.0005538235099743024
                                     0.001588211190322838
                                    -0.010594910698770119
                                    -0.0003723971024834356]

    @test gnr_fes_res["E"] == 1.0360429884274056

    @test Set(gnr_fes_res["polynom_series"]) == Set([:k, :k_k; :k_k_k; :k_k_i; :k_i; :k_i_i; :i; :i_i; :i_i_i])

    @test Set(gnr_fes_res["all_inputs"]) == Set([:k; :i])

    @test gnr_fes_res["fixed_inputs"] == [:k]

    @test gnr_fes_res["flexible_input"] == :i

    @test gnr_fes_res["share_degree"] == 3

    @test size(gnr_fes_res["fes_optim_estimates"]) == (10, 4)

    @test typeof(gnr_fes_res["fes_optim_estimates"]) == DataFrame

    @test names(gnr_fes_res["fes_optim_estimates"]) == ["Variable"; "γ"; "γ_dash"; "γ_flex"]

    @test Optim.converged(gnr_fes_res["fes_optim_results"])

    @test gnr_ses_res["α"] == [-0.024846888000800222; 0.3882501166382591; 0.002466151018639108]

    @test gnr_ses_res["δ"] == [0.16802673268297028; 0.7690654898797339; 0.06544391037271638; -0.03980526457212217;;]
    
    @test Optim.converged(gnr_ses_res["gmm_optim_results"]) 

    @test Set(gnr_ses_res["polynom_fixed"]) == Set([:k; :k_k; :k_k_k])

    @test gnr_ses_res["fixed_inputs"] == [:k]

    @test gnr_ses_res["flexible_input"] == :i

    @test gnr_ses_res["int_const_series_degree"] == 3

    @test gnr_ses_res["lm_tfp_degree"] == 3

    @test Set(gnr_ses_res["all_inputs"]) == Set([:k; :i])
end