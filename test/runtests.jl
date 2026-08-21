# Load required dependencies (use 'using Pkg' and Pkg.add("DataFrames"), Pkg.add("GNRProdEst") or Pkg.add("CSV") if you do not have any of these packages in your environment.)
# The actual tests live in the @testitem blocks in this folder (see estimation_tests.jl).
# They are discovered by the VS Code test explorer directly, and by Pkg.test through the runner below.
using TestItemRunner

@run_package_tests
