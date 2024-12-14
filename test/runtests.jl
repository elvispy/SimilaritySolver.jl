using SimilaritySolver
using Test
using Symbolics

@testset "Dependence tests" begin
    @variables x y η(x, y) f(η)
    @test SimilaritySolver.is_output(η, [x, y]) == true 
    @test SimilaritySolver.is_output(f, [x, y]) == true
    @test SimilaritySolver.is_output(f, [η, f]) == true 
    @test SimilaritySolver.is_output(x, [x, y]) == false 
    @test SimilaritySolver.is_output(y, [x, f]) == false 


    @test SimilaritySolver.is_input(η, [x, y]) == false 
    @test SimilaritySolver.is_input(f, [x, y]) == false
    @test SimilaritySolver.is_input(f, [η, f]) == false 
    @test SimilaritySolver.is_input(x, [x, y]) == false 
    @test SimilaritySolver.is_input(y, [x, f]) == true 
    @test SimilaritySolver.is_input(η, [x, f]) == true 
    @test SimilaritySolver.is_input(η, [x, f]) == true 

end

@testset "Identifying inputs/outputs" begin
    @variables ν x y u(x, y)
    eqn = Differential(x)(u) - ν * Differential(y)(Differential(y)(u))
    vars = Num.(Symbolics.get_variables(eqn));
    ins, outs = SimilaritySolver.decomposeVars(vars ∪ [x, y]);
    @test isequal(ins)([x, y])
    @test isequal(outs)([u])
end

# @testset "SimilaritySolver.jl" begin
#     @test 
# end
