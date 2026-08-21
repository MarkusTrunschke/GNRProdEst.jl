# Shared setup: estimating both stages is expensive, so it is done once here and the
# results are reused by every @testitem that declares setup=[ReplicationResults].
@testmodule ReplicationResults begin

    using GNRProdEst, DataFrames, CSV, Optim

    # Read in replication data. @__DIR__ is this file's directory, so the data is found no
    # matter what the working directory is.
    rep_data = CSV.read(joinpath(@__DIR__, "GNR_data_500.csv"), DataFrame)

    # The reference values asserted below were produced with these GMM tolerances, so this
    # dict has to stay as it is. Printing is disabled: no test inspects the printed tables and
    # it only slows the suite down.
    opts = Dict("ses_optimizer_options" => Optim.Options(f_tol = 1e-12,
                                                         x_tol = 1e-12,
                                                         g_tol = 1e-12),
                "fes_print_starting_values" => false,
                "ses_print_starting_values" => false,
                "fes_print_results" => false,
                "ses_print_results" => false,
                "print_results" => false
                )

    # Run both estimation stages at the same time. boot_reps is kept small on purpose: no test
    # below asserts on standard errors, and every repetition re-estimates the whole model.
    gnr_fes_res, gnr_ses_res = GNRProdEst.gnrprodest(data = rep_data,
                                                     output = :yg,
                                                     flexible_input = :i,
                                                     fixed_inputs = :k,
                                                     ln_share_flex_y = :si,
                                                     id = :id,
                                                     time = :time,
                                                     boot_reps = 2,
                                                     opts = opts);
end
