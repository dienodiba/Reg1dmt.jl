using Reg1dmt
using Test

@testset "Reg1dmt.jl" begin
    @test Reg1dmt.Inversion("input_data.dat","input_setting.dat","output.dat") == nothing
end