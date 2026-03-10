# SimilaritySolver.jl

*Automated Similarity Reduction for PDEs in Julia*

`SimilaritySolver.jl` is a high-performance symbolic tool for reducing Partial Differential Equations (PDEs) to Ordinary Differential Equations (ODEs) using exact scaling (dilation) symmetry.

## Overview

Finding similarity solutions traditionally requires manual identification of scaling invariances—a tedious process of "guessing" exponents that balance terms in a PDE. `SimilaritySolver.jl` automates this process through:

1.  **Dilation Symmetry Detection**: Uses exact rational linear algebra to find the null space of the scaling-invariance system.
2.  **PDE→ODE Reduction**: Automatically derives the similarity variable $\eta$ and transforms the PDE into a reduced ODE system.
3.  **Symbolic Integration**: Built on `Symbolics.jl` for robust mathematical manipulation.

## Key Features

- **Exact Arithmetic**: Computes scaling exponents as exact rationals ($1/2$, $2/3$, etc.), avoiding numerical drift.
- **System Support**: Handles coupled systems of PDEs (e.g., fluid-surfactant interactions).
- **Parsimony Ordering**: Returns the simplest similarity variables first (e.g., $\eta = x/\sqrt{t}$).
- **String & Symbolic APIs**: Use natural string notation or full `Symbolics.Num` expressions.

## Installation

```julia
using Pkg
Pkg.add("SimilaritySolver")
```

## Quick Example

```julia
using SimilaritySolver, Symbolics

# Heat equation: u_t = u_xx
@variables x t u(..)
Dt = Differential(t); Dx = Differential(x)
heat = Dt(u(x,t)) - Dx(Dx(u(x,t)))

results = find_ode_dilation(heat; indep_vars=[x,t], dep_vars=[u(x,t)])

# View the similarity variable
println(results[1]["similarity_variable"])  # x*t^(-1//2)
```

## Community

`SimilaritySolver.jl` is a community-driven project. We welcome bug reports, feature requests, and pull requests. Please see the [Contributing](contributing.md) guide for more details.
