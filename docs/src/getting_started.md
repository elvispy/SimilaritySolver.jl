# Getting Started

Similarity solutions are a class of invariant solutions for Partial Differential Equations (PDEs). They allow us to transform a PDE with $n$ independent variables into an Ordinary Differential Equation (ODE) or a PDE with $n-1$ variables.

## Installation

`SimilaritySolver.jl` can be installed using the Julia package manager:

```julia
using Pkg
Pkg.add("SimilaritySolver")
```

Alternatively, to install the latest development version:

```julia
Pkg.add(url="https://github.com/elvispy/SimilaritySolver.jl.git")
```

## Verifying the Installation

After installation, it is recommended to run the test suite to ensure everything is working correctly on your system:

```julia
Pkg.test("SimilaritySolver")
```

All 280+ tests should pass. If you encounter any failures, please see the [Contributing](contributing.md) guide to report them.

## Basic Workflow

1.  **Define your PDE**: You can use either Julia `Symbolics.jl` expressions or a natural string notation.
2.  **Call the solver**: Use `find_ode_dilation` (for symbolic input) or `find_similarity_v2` (for string input).
3.  **Inspect results**: The solver returns a dictionary containing the similarity variable $\eta$, the reduced ODE, and the dependent-variable ansatz.

## Example 1: Heat Equation (Symbolic)

```julia
using SimilaritySolver, Symbolics

# 1. Define symbolic variables and operators
@variables x t u(..)
Dt = Differential(t); Dx = Differential(x)

# 2. Set up the heat equation: u_t - u_xx = 0
heat = Dt(u(x,t)) - Dx(Dx(u(x,t)))

# 3. Solve for similarity transformations
results = find_ode_dilation(heat; indep_vars=[x,t], dep_vars=[u(x,t)])

# 4. Extract first candidate
res = results[1]
println("Similarity Variable η: ", res["similarity_variable"])
println("Reduced ODE: ", res["PDE_similarity"])
```

## Example 2: KdV Equation (String API)

The string API is convenient for quick experiments or interactive use.

```julia
using SimilaritySolver

# Korteweg-de Vries (KdV): u_t + 6u·u_x + u_xxx = 0
pde = "du/dt + 6*u*du/dx + d3u/d3x = 0"
bcs = "u(x=Inf, t) = 0"

result = find_similarity_v2(pde, bcs)

println("η = ", result["similarity_variable"])
println("Reduced ODE: ", result["PDE_similarity"])
```

## Coupled Systems

`SimilaritySolver.jl` can also handle systems of PDEs. Simply pass a vector of equations to `find_ode_dilation`.

```julia
# IM-regime from EDS 2025
pdes = [mass, momentum, surfactant]
results = find_ode_dilation(pdes; indep_vars=[r,t], dep_vars=[u(r,t), h(r,t), Γ(r,t)])
```
