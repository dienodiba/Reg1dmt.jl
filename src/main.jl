# Dieno Diba 2022
# REG1DMT
# Main program

using DelimitedFiles
using Statistics
using SparseArrays
using LinearAlgebra
using Printf
using Dates

function main()

    filedata = readline()
    filestg = readline()
    fileout = readline()

    Inv1DMT(filedata,filestg,fileout)
    
end

include("Inv1DMT.jl")
include("Fwd1DMT.jl")

main()
