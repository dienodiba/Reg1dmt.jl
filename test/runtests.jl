using REG1DMT
using Test

@testset "REG1DMT.jl" begin
    @test REG1DMT.Inversion("input_data.dat","input_setting.dat","output.dat") == nothing
end