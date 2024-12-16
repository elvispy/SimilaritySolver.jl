using SimilaritySolver
using Test
using Symbolics
using Combinatorics

iszerosym(sym) = (simplify(sym; expand=true) != 0) === false

@testset "parse_boundary_conditions" begin
    
    #@variables x t u(t, x)# η(x, y) f(η)
    boundary_conditions = "w(r=Inf, t) = -Inf; w(r, t=0) = 1.0"
    boundary_list = String.(split(boundary_conditions, ";"))

    boundary_conditions = SimilaritySolver.parse_boundary_condition.(boundary_list)
    #println(boundary_conditions[1])
    input_vars = convert(Vector{Num}, union(map(var -> var[2], boundary_conditions)...));
    restrictions = map(var -> var[1], boundary_conditions);
    output_vars = convert(Vector{Num}, union(map(var -> var["function"], restrictions)...));
    ss = @variables r, t w(r, t)
    @test isequal(Set(input_vars))(Set(ss))
    @test isinf(boundary_conditions[1][1]["restriction"][r])
    @test isnothing(boundary_conditions[1][1]["restriction"][t])
    @test isnothing(boundary_conditions[2][1]["restriction"][r])
    @test iszero(boundary_conditions[2][1]["restriction"][t])
    @test isinf(boundary_conditions[1][1]["value"]) && boundary_conditions[1][1]["value"] < 0 
    @test boundary_conditions[2][1]["value"] ≈ 1.0

end

@testset "is_output/is_input" begin
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

@testset " decomposeVars    " begin
    @variables ν x y u(x, y) v(x, y) β
    eqn = Differential(x)(u) - ν * Differential(y)(Differential(y)(u))
    vars = Num.(Symbolics.get_variables(eqn));
    ins, outs = SimilaritySolver.decomposeVars(vars ∪ [x, y]);
    #@test isequal(ins)([x, y])
    @test isequal(outs)([u])

    vals = [ν, x, y, u, v, β]
    belongs(value, list) = any(isequal(value), list);
    @testset "Testing outputs in $subset" for subset = filter(subset -> length(subset) >= 3, collect(powerset(vals)))
        ins, outs = SimilaritySolver.decomposeVars(subset)
        if belongs(u, subset)
            @test belongs(u, outs) == any(map(val -> belongs(val, subset), [x, y]))
            @test belongs(u, ins)  == false #du/dv = 0
        end
        if belongs(v, subset)
            @test belongs(v, outs) == any(map(val -> belongs(val, subset), [x, y]))
            @test belongs(v, ins)  == false

        end
        if belongs(β, subset)
            @test belongs(β, outs) == false
            @test belongs(β, ins)  == false
        end
    end
end

@testset "find_ode/heatEqn  " begin
    @variables w k u(w, k)# η(x, y) f(η)
    symbolicPDE = Differential(w)(u) - Differential(k)(Differential(k)(u))
    vars = [w, k, u];
    result = SimilaritySolver.find_ode(symbolicPDE; vars=vars, log=false);
    @test result !== nothing
    # @test isequal(Set(values(result["substitutions"])))(Set([0., -0.5]))
    # η, f = Symbolics.get_variables(result["PDE"]; sort = true)

    # if !iszerosym(Symbolics.expand_derivatives(Differential(f)(η))); η, f = f, η; end
    # @test iszerosym(Symbolics.coeff(result["PDE"], Differential(η)(f))+.5* η)
    # @test iszerosym(Symbolics.coeff(result["PDE"], Differential(η)(Differential(η)(f)))+ 1.0)
    # expr = result["PDE"] - (-0.5*η*Differential(η)(f) - Differential(η)(Differential(η)(f)));
    
    # @test iszerosym(expr)

end

 @testset "find_ode/kdV      " begin
     # Source: https://www.ucl.ac.uk/~ucahhwi/LTCC/sectionB-similarity.pdf
     @variables x t u(t, x)# η(x, y) f(η)
     Dt = Differential(t); Dx = Differential(x)
     symbolicPDE = Dt(u) + 6*u*Dx(u) + Dx(Dx(Dx(u)))
     vars = [t, x, u];
     result = SimilaritySolver.find_ode(symbolicPDE; vars=vars, log=false)
     #println(result)
     @test result !== nothing
     #@test isequal(Set(values(result["substitutions"])))(Set([1.0, -0.5]))
     #η, f = Symbolics.get_variables(result["PDE"]; sort = true)
     #if !iszerosym(Symbolics.expand_derivatives(Differential(f)(η))); η, f = f, η; end
     #@test iszerosym(Symbolics.coeff(result["PDE"], Differential(η)(f))+2.0+ .5* η^2)
     #@test iszerosym(Symbolics.coeff(result["PDE"], Differential(η)(Differential(η)(f)))+ η^2)

     #expr = result["PDE"] - (-2.0*Differential(η)(f) - η*Differential(η)(Differential(η)(f)) - 0.5*(η^2)*Differential(η)(f));

     #@test iszerosym(expr)

 end

