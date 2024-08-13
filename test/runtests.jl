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
                                                            boot_reps = 10,
                                                            opts = opts
                                                    );

    @test [0.6523; -0.0012;  0.0046; -0.0015;  0.001; -4.6e-5; -0.0005;  0.0012; -0.0004; -0.022] < gnr_fes_res["γ"] <  [0.6524; -0.0001;  0.0048; -0.0014;  0.0012; -4.4e-5; -0.0004;  0.0013; -0.0002; -0.0211]

    @test [ 0.675; -0.002;  0.004; -0.002;  0.001; -4.7e-5; -0.0006;  0.001; -0.002; -0.03] < gnr_fes_res["γ_dash"] < [ 0.676; -0.001;  0.005; -0.001;  0.002; -4.5e-5; -0.0004;  0.002; -0.0001; -0.02]

    @test [ 0.65; -0.0001;  0.001; -0.002;  0.0001; -1.9e-5; -0.0009;  0.01; -0.001; -0.02] <gnr_fes_res["γ_flex"] < [ 0.66; -0.001;  0.002; -0.001;  0.001; -1.0e-5; -0.0001;  0.001; -0.0001; -0.005]

    @test 1.035 < gnr_fes_res["E"] < 1.037

    @test Set(gnr_fes_res["polynom_series"]) == Set(Symbol.(["k"; "k⋅k"; "k⋅k⋅k"; "k⋅k⋅i"; "k⋅i"; "k⋅i⋅i"; "i"; "i⋅i"; "i⋅i⋅i"]))

    @test Set(gnr_fes_res["all_inputs"]) == Set([:k; :i])

    @test gnr_fes_res["fixed_inputs"] == [:k]

    @test gnr_fes_res["flexible_input"] == :i

    @test gnr_fes_res["share_degree"] == 3

    @test size(gnr_fes_res["fes_optim_estimates"]) == (10, 4)

    @test typeof(gnr_fes_res["fes_optim_estimates"]) == DataFrame

    @test names(gnr_fes_res["fes_optim_estimates"]) == ["Variable"; "γ"; "γ_dash"; "γ_flex"]

    @test Optim.converged(gnr_fes_res["fes_optim_results"])

    @test [0.38, -0.02, 0.0001] < gnr_ses_res["α"] < [0.4, -0.03, 0.004]

    @test [0.16; 0.76; 0.06; -0.05] < vec(gnr_ses_res["δ"]) < [0.175; 0.775; 0.07; -0.05]
    
    @test Optim.converged(gnr_ses_res["gmm_optim_results"]) 

    @test Set(gnr_ses_res["polynom_fixed"]) == Set(Symbol.(["k"; "k⋅k"; "k⋅k⋅k"]))

    @test gnr_ses_res["fixed_inputs"] == [:k]

    @test gnr_ses_res["flexible_input"] == :i

    @test gnr_ses_res["int_const_series_degree"] == 3

    @test gnr_ses_res["lm_tfp_degree"] == 3

    @test Set(gnr_ses_res["all_inputs"]) == Set([:k; :i])
end
