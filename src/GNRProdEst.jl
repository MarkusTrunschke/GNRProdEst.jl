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
    using DataFrames, Revise, NonlinearSolve, GLM, LsqFit, PrettyTables, Statistics

    ## Include other files defining functions
    include("data_prep.jl") # Set up inputs to be in the correct for for the internals to work with
    include("estimation_fnc.jl") # Function to estimate the production function

    export Nothing

end