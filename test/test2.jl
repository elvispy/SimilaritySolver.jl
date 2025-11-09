using SimilaritySolver

# ν and U∞ are treated as symbols/parameters
pde = "dψ/dy * d2ψ/dxdy - dψ/dx * d2ψ/d2y - ν * d3ψ/d3y = 0"

# boundary conditions: ψ(x,0)=0, ψ_y(x,0)=0, ψ_y(x,∞)=U∞
bcs = "ψ(x, y=0) = 0; dψ/dy(x, y=0) = 0; dψ/dy(x, y=Inf) = U∞"

result = find_similarity(pde, bcs; parameters=["ν", "U∞"])
println(result)            # shows inputs/outputs, similarity variable η, reduced ODE, and transformed BCs



# Test the function with an example
#parsed_conditions = find_similarity("dw/dt = w * d2w/d2r", "w(r=Inf, t) = 0; w(r, t=0) = U0")
#parsed_conditions = find_similarity("du/dt + 6 * u * du/dx + d3u/d3x = 0", "u(x=Inf, t) = 0")
#println(parsed_conditions)

# @variables t r u(r, t)# η(x, y) f(η)
# symbolicPDE = Differential(t)(u) - Differential(r)(Differential(r)(u))
# vars = [t, r, u];
# result = find_ode(symbolicPDE; vars=vars)
# @variables x t u(t, x)# η(x, y) f(η)
# Dt = Differential(t); Dx = Differential(x)
# symbolicPDE = Dt(u) + 6*u*Dx(u) + Dx(Dx(Dx(u)))
# vars = [t, x, u];
# result = find_ode(symbolicPDE; vars=vars)
