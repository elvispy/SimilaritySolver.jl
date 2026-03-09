using SimilaritySolver, LinearAlgebra
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
    println(boundary_conditions)
    @test isequal(Set(input_vars))(Set(ss))
    @test isinf(boundary_conditions[1][1]["restriction"][r])
    @test !haskey(boundary_conditions[1][1]["restriction"], t)
    @test !haskey(boundary_conditions[2][1]["restriction"], r)
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

@testset "find_ode_dilation/heatEqn" begin
    @variables x t u(..)
    Dt = Differential(t); Dx = Differential(x)
    heat = Dt(u(x,t)) - Dx(Dx(u(x,t)))
    results = find_ode_dilation(heat; indep_vars=[x,t], dep_vars=[u(x,t)])
    @test !isempty(results)
    # similarity variable should be η = x / t^(1/2)
    η_expr = results[1]["similarity_variable"]
    @test occursin("t", string(η_expr))  # t appears in η
    @test occursin("x", string(η_expr))  # x appears in η
    ode_str = string(results[1]["PDE"])
    # ODE must contain f'' (2nd derivative of f_dil)
    @test occursin("Differential(η_bare)(Differential(η_bare)(f_dil", ode_str)
    # alpha_vec and gamma are scaling exponents (rational)
    @test results[1]["alpha_vec"] isa Vector{<:Rational}
    @test results[1]["gamma"] isa Rational
end

@testset "find_ode_dilation/KdV" begin
    @variables x t u(..)
    Dt = Differential(t); Dx = Differential(x)
    kdv = Dt(u(x,t)) + 6*u(x,t)*Dx(u(x,t)) + Dx(Dx(Dx(u(x,t))))
    results = find_ode_dilation(kdv; indep_vars=[x,t], dep_vars=[u(x,t)])
    @test !isempty(results)
    # similarity variable: η = x / t^(1/3)
    η_str = string(results[1]["similarity_variable"])
    @test occursin("1//3", η_str) || occursin("(1//3)", η_str)
    # scaling exponent: γ = -2
    @test results[1]["gamma"] == -2//1
    ode_str = string(results[1]["PDE"])
    # ODE must contain f''' (3rd derivative of f_dil)
    @test occursin("Differential(η_bare)(Differential(η_bare)(Differential(η_bare)(f_dil", ode_str)
end

@testset "find_ode_dilation/3vars" begin
    # 2D heat equation: u_t = u_xx + u_yy  (3 independent variables: x, y, t)
    # Dilation symmetry: (t,x,y) ~ (λ²t, λx, λy), u invariant.
    # Similarity variables: η = x/(√t·y), η = y/(√t·x), or η = t/(x²y²).
    @variables x y t u(..)
    Dt = Differential(t); Dx = Differential(x); Dy = Differential(y)
    heat2d = Dt(u(x,y,t)) - Dx(Dx(u(x,y,t))) - Dy(Dy(u(x,y,t)))
    results = find_ode_dilation(heat2d; indep_vars=[x,y,t], dep_vars=[u(x,y,t)], verbose=true)
    # Architecture test: pipeline handles 3 independent variables without error
    @test !isempty(results)
    # alpha_vec has 3 entries (one per independent variable)
    @test length(results[1]["alpha_vec"]) == 3
    # At least one similarity variable involves all three variables
    η_strings = [string(r["similarity_variable"]) for r in results]
    @test any(s -> occursin("x", s) && occursin("y", s) && occursin("t", s), η_strings)
    # All ODEs are free of the original independent variables
    for r in results
        ode_str = string(r["PDE"])
        @test !occursin("u(x, y, t)", ode_str)
    end
end

@testset "find_ode_dilation/Burgers" begin
    # u_t + u*u_x = ν*u_xx  (ν treated as fixed parameter)
    # Dilation: x~λx, t~λ²t, u~λ⁻¹u  →  η = x/√t, γ = -1
    # ODE: -F/2 - η·F'/2 + F·F' = ν·F''  (Hopf-Cole structure)
    @variables x t u(..) ν
    Dt = Differential(t); Dx = Differential(x)
    burgers = Dt(u(x,t)) + u(x,t)*Dx(u(x,t)) - ν*Dx(Dx(u(x,t)))
    results = find_ode_dilation(burgers; indep_vars=[x,t], dep_vars=[u(x,t)])
    @test !isempty(results)
    # scaling exponent γ = -1
    @test results[1]["gamma"] == -1//1
    # similarity variable contains x and t (η = x/√t)
    η_str = string(results[1]["similarity_variable"])
    @test occursin("x", η_str) && occursin("t", η_str)
    # ODE must contain f'' (viscous diffusion term)
    ode_str = string(results[1]["PDE"])
    @test occursin("Differential(η_bare)(Differential(η_bare)(f_dil", ode_str)
end


@testset "find_ode_dilation/EDS2025_IM" begin
    # Eshima, Deike & Stone (J. Fluid Mech. 2025) — inertio-Marangoni (IM) regime
    # Coupled PDEs for (u, h, Γ) in radial coordinates (eqs. 1.3a–c):
    #   ∂h/∂t + (1/r)∂(r·u·h)/∂r = 0                   (mass)
    #   ∂u/∂t + u·∂u/∂r = -(2/h)·∂Γ/∂r                 (momentum, σ = -Γ)
    #   ∂Γ/∂t + (1/r)∂(r·u·Γ)/∂r = 0                   (surfactant)
    #
    # Dilation: (r,t,u,h,Γ) ~ (λ³,λ⁴,λ⁻¹,λ⁰,λ⁻²), i.e. null vec [3,4,-1,0,-2]
    #   → η = r/t^(3/4),   u ~ t^(-1/4),  h ~ t^0,  Γ ~ t^(-1/2)  (EDS Table 1)
    @variables r t u(..) h(..) Γ(..)
    Dr = Differential(r); Dt = Differential(t)

    mass        = Dt(h(r,t)) + (1/r)*Dr(r*u(r,t)*h(r,t))
    momentum    = Dt(u(r,t)) + u(r,t)*Dr(u(r,t)) - (2/h(r,t))*Dr(Γ(r,t))
    surfactant  = Dt(Γ(r,t)) + (1/r)*Dr(r*u(r,t)*Γ(r,t))
    pdes = [mass, momentum, surfactant]

    # 1. Invariance system: 5 scaling params [a_r, a_t, c_u, c_h, c_Γ], rank 2
    A, param_names = SimilaritySolver.build_invariance_system(
        pdes, [r,t], [u(r,t), h(r,t), Γ(r,t)])
    @test size(A, 2) == 5          # 2 indep + 3 dep params
    @test rank(A) == 2             # two independent constraints → 3D null space

    # 2. EDS null vector [a_r=3, a_t=4, c_u=-1, c_h=0, c_Γ=-2] satisfies invariance
    r_idx  = findfirst(==(:a_r), param_names)
    t_idx  = findfirst(==(:a_t), param_names)
    u_idx  = findfirst(==(:c_u), param_names)
    h_idx  = findfirst(==(:c_h), param_names)
    Γ_idx  = findfirst(==(:c_Γ), param_names)
    eds_vec = zeros(5)
    eds_vec[r_idx]=3; eds_vec[t_idx]=4; eds_vec[u_idx]=-1; eds_vec[h_idx]=0; eds_vec[Γ_idx]=-2
    @test all(abs.(A * eds_vec) .< 1e-10)

    # 3. Direct reduction with the known EDS null vector.
    #    (full pipeline search is O(N³) combinations and takes >200s; tested separately.)
    eds_null = Rational{Int}[3, 4, -1, 0, -2]
    cands = SimilaritySolver.build_similarity_candidates(
        eds_null, param_names, [r,t], [u(r,t), h(r,t), Γ(r,t)])
    @test length(cands) == 2   # pivot=r and pivot=t

    # Pick pivot=t candidate: γ_u=-1/4, γ_h=0, γ_Γ=-1/2 (EDS Table 1)
    t_cand = cands[findfirst(c -> isequal(c.pivot, t), cands)]
    @test t_cand.gamma_vals == [-1//4, 0//1, -1//2]

    # 4. reduce_to_ode produces a system of 3 ODEs free of r and t
    ode_sys = SimilaritySolver.reduce_to_ode(
        pdes, [r,t], [u(r,t), h(r,t), Γ(r,t)],
        t_cand.η_expr, t_cand.pivot, t_cand.pivot_idx,
        t_cand.gamma_vals, t_cand.alpha_norm)
    @test ode_sys !== nothing
    @test ode_sys isa AbstractVector
    @test length(ode_sys) == 3    # one ODE per PDE in the system
    for ode in ode_sys
        ode_str = string(ode)
        @test !occursin("u(r, t)", ode_str)
        @test !occursin("h(r, t)", ode_str)
        @test !occursin("Γ(r, t)", ode_str)
    end
end

@testset "find_ode_dilation/nonlinear_diffusion" begin
    # Nonlinear diffusion (n=2): u_t = (u² u_x)_x  — Barenblatt-Pattle equation
    # Known dilation symmetry: (x,t,u) ~ (λx, λ⁴t, λ⁻¹u)
    # i.e., scaling vector (a_x=1, a_t=4, a_u=-1) → similarity variable η = x/t^(1/4), γ = -1/4
    @variables x t u(..)
    Dt = Differential(t); Dx = Differential(x)
    pde = Dt(u(x,t)) - Dx(u(x,t)^2 * Dx(u(x,t)))

    # 1. Invariance system: one constraint (rank 1), three scaling params
    A, param_names = SimilaritySolver.build_invariance_system(pde, [x,t], [u(x,t)])
    @test size(A, 2) == 3          # params: a_x, a_t, a_u
    @test rank(A) == 1             # one independent constraint → 2D null space

    # 2. Barenblatt null vector [a_x=1, a_t=4, a_u=-1] satisfies the invariance system.
    #    The exponent a_t/a_x = 4 = n+2 encodes the known 1/(n+2) = 1/4 similarity exponent.
    x_idx = findfirst(==(:a_x), param_names)
    t_idx = findfirst(==(:a_t), param_names)
    u_idx = findfirst(s -> s != :a_x && s != :a_t, param_names)
    barenblatt = zeros(3); barenblatt[x_idx]=1; barenblatt[t_idx]=4; barenblatt[u_idx]=-1
    @test all(abs.(A * barenblatt) .< 1e-10)

    # 3. Pipeline finds the Barenblatt solution via null-space combination search (Phase 2.5).
    #    The combination 2v₁ - v₂ of the two basis vectors yields α=(1,4,-1).
    results = find_ode_dilation(pde; indep_vars=[x,t], dep_vars=[u(x,t)], verbose=true)
    @test !isempty(results)
    @test all(r -> r["alpha_vec"] isa Vector{<:Rational}, results)
    # Barenblatt: γ = -1/4, η = x / t^(1/4)
    @test any(r -> r["gamma"] == -1//4, results)
    η_strs = [string(r["similarity_variable"]) for r in results]
    @test any(s -> occursin("1//4", s) || occursin("(1//4)", s), η_strs)
end

