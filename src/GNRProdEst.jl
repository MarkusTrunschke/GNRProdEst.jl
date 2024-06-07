#####################################################
############### THE GNRProdEst PACKAGE ##############
# Author: Markus Trunschke                          #
# Date of first version: 18.04.2024                 #
# Contact: markus.trunschke@googlemail.com          #
#####################################################
# Notes: This package implements the estimation of  #
#        a production function using the method de- #
#        scribed in Gandhi, Navarro, Rivers (2020). #
#####################################################

module GNRProdEst
    # Load dependencies
    using DataFrames, Revise, NonlinearSolve, PrettyTables, Statistics, ShiftedArrays, LinearAlgebra, Optim, BenchmarkTools, NaNMath
    
    ## Set a number of parameters
    # Set BLAS and LAPACK function number of threads to 1 because it interferes with bootstrap loop parallelization. # Why? Well, what do I know? I just noticed that if I run a lot of stuff in parallel (EVEN THOUGH IT USES A LOT OF BLAS OPERATIONS!),
    BLAS.set_num_threads(1)

    ## Include other files defining functions
    include("data_prep.jl") # Set up inputs to be in the correct for for the internals to work with
    include("estimation_fnc.jl") # Function to estimate the production function
    include("aux_fnc.jl") # All auxiliary function definitions

    # Export user functions
    export gnrprodest, gnrprodest!, gnrfirststage, gnrfirststage!, gnrsecondstage, gnrseconstage!

end