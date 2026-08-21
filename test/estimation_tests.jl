@testitem "First stage: coefficient estimates" setup=[ReplicationResults] begin

    gnr_fes_res = ReplicationResults.gnr_fes_res

    @test round.(gnr_fes_res["γ"], digits = 5) == [0.65239, -0.00112, 0.00476, -0.00146, 0.00111, -5.0e-5, -0.0005, 0.00123, -0.00033, -0.02119]

    @test round.(gnr_fes_res["γ_dash"], digits = 5) == [0.6759, -0.00116, 0.00494, -0.00152, 0.00115, -5.0e-5, -0.00052, 0.00128, -0.00035, -0.02195]

    @test round.(gnr_fes_res["γ_flex"], digits = 5) == [0.65239, -0.00037, 0.00159, -0.00146, 0.00055, -1.0e-5, -0.0005, 0.00062, -0.00033, -0.01059]

    @test gnr_fes_res["E"] ≈ 1.03604299 rtol = 1e-8
end

@testitem "First stage: series and inputs" setup=[ReplicationResults] begin

    gnr_fes_res = ReplicationResults.gnr_fes_res

    @test Set(gnr_fes_res["polynom_series"]) == Set(Symbol.(["k"; "k⋅k"; "k⋅k⋅k"; "k⋅k⋅i"; "k⋅i"; "k⋅i⋅i"; "i"; "i⋅i"; "i⋅i⋅i"]))

    @test Set(gnr_fes_res["all_inputs"]) == Set([:k; :i])

    @test gnr_fes_res["fixed_inputs"] == [:k]

    @test gnr_fes_res["flexible_input"] == :i

    @test gnr_fes_res["share_degree"] == 3
end

@testitem "First stage: optimization output" setup=[ReplicationResults] begin

    using DataFrames, Optim

    gnr_fes_res = ReplicationResults.gnr_fes_res

    @test size(gnr_fes_res["fes_optim_estimates"]) == (10, 4)

    @test typeof(gnr_fes_res["fes_optim_estimates"]) == DataFrame

    @test names(gnr_fes_res["fes_optim_estimates"]) == ["Variable"; "γ"; "γ_dash"; "γ_flex"]

    @test Optim.converged(gnr_fes_res["fes_optim_results"])
end

@testitem "Second stage: coefficient estimates" setup=[ReplicationResults] begin

    using Optim

    gnr_ses_res = ReplicationResults.gnr_ses_res

    @test round.(gnr_ses_res["α"], digits = 5) == [0.38825, -0.02485, 0.00247]

    # 4 digits: δ[4] differs in the 5th digit between a plain run and one under
    # --check-bounds=yes (which is what Pkg.test uses)
    @test round.(vec(gnr_ses_res["δ"]), digits = 4) == [0.168, 0.7691, 0.0654, -0.0398]

    @test Optim.converged(gnr_ses_res["gmm_optim_results"])
end

@testitem "Second stage: series and inputs" setup=[ReplicationResults] begin

    gnr_ses_res = ReplicationResults.gnr_ses_res

    @test Set(gnr_ses_res["polynom_fixed"]) == Set(Symbol.(["k"; "k⋅k"; "k⋅k⋅k"]))

    @test gnr_ses_res["fixed_inputs"] == [:k]

    @test gnr_ses_res["flexible_input"] == :i

    @test gnr_ses_res["int_const_series_degree"] == 3

    @test gnr_ses_res["lm_tfp_degree"] == 3

    @test Set(gnr_ses_res["all_inputs"]) == Set([:k; :i])
end
