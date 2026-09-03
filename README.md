[![Build Status](https://github.com/elvispy/SimilaritySolver.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/elvispy/SimilaritySolver.jl/actions/workflows/CI.yml?query=branch%3Amain)

# Symbolic Similarity PDE Solver (Julia)

This Julia package provides a set of tools for symbolic manipulation and similarity transformation of partial differential equations (PDEs). It is built on top of the `Symbolics.jl` and `SymbolicUtils.jl` ecosystems and supports automatic classification of variables, similarity analysis, boundary condition parsing, and ODE reduction.

---

## Motivation

Similarity solutions are a powerful technique for reducing partial differential equations (PDEs) to ordinary differential equations (ODEs). This simplification often makes complex analytical or numerical problems more tractable. However, identifying appropriate transformations and substitutions is nontrivial.

This package automates that process by:
- Classifying symbolic variables as inputs or outputs
- Substituting similarity variables
- Applying heuristics to guess transformation powers
- Verifying and simplifying the reduced equation

This tool aims to accelerate mathematical modeling in physics, engineering, and applied mathematics.

---

## Features

- **Automatic decomposition** of symbolic variables into inputs and outputs.
- **Similarity transformation detection** for 2D PDEs.
- **Boundary condition parsing** and transformation under similarity variables.
- **Power guessing heuristic** for identifying simplifying substitutions.
- **Parsing PDEs from strings** using symbolic differential notation like `d2x/dy`.

---

## 📦 Installation

Make sure you have Julia installed. Then clone this repository and instantiate dependencies as such:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Dependencies include:
- `Symbolics`
- `SymbolicUtils`
- `Logging`

---

## 🔧 Usage Example

### Find Similarity Solution

```julia
using Symbolics

result = find_similarity("du/dt + 6 * u * du/dx + d3u/d3x = 0", "u(x=Inf, t) = 0")
println(result)
```

### Blasius Boundary Layer

The classical flat-plate boundary layer satisfies

```text
ψ_y ψ_{xy} - ψ_x ψ_{yy} - ν ψ_{yyy} = 0
```

with boundary conditions `ψ(x,0)=0`, `∂ψ/∂y(x,0)=0`, and `∂ψ/∂y(x,∞)=U∞`. You can reproduce
the Blasius similarity reduction with:

```julia
using SimilaritySolver

pde = "dψ/dy * d2ψ/dxdy - dψ/dx * d2ψ/d2y - ν * d3ψ/d3y = 0"
bcs = "ψ(x, y=0) = 0; dψ/dy(x, y=0) = 0; dψ/dy(x, y=Inf) = U∞"

result = find_similarity(pde, bcs; parameters=["ν", "U∞"])
println(result["similarity_variable"])  # η(x, y) => y * x^m guess that succeeded
println(result["output_similarity"])    # ψ(x, y) => x^n f(η)
println(result["PDE_similarity"])       # reduced ODE (Blasius: f''' + 0.5 f f'' = 0)
```

`result["PDE_similarity"]` contains the familiar Blasius ODE (up to scaling by the supplied parameters), making it easy to hand the reduced system to an ODE solver.

### Parse Boundary Conditions

```julia
bc = "f(x=0, t) = 1.0"
parsed = parse_boundary_condition(bc)
println(parsed[1]["function"])
println(parsed[1]["restriction"])
println(parsed[1]["value"])
```

### Parse a PDE String

```julia
parse_pde("d2x/dy = x * y", [:x, :y], [:u])
```

---

## 🔤 Notation and Syntax

The core function `find_ode` and its wrapper `find_similarity` expect PDEs and boundary conditions in a symbolic or string-based format. Here’s how the notation works:

### PDE Expression Syntax
- Use `d<order><variable>` for partial derivatives.
- Use `/` to indicate the variable of differentiation.
- Use `=` to denote equality (internally rewritten as `- RHS`).

Example:
```julia
"du/dt + 6 * u * du/dx + d3u/d3x = 0"
```
will be interpreted as:
$\frac{\partial u}{\partial t} + 6u \frac{\partial u}{\partial x} + \frac{\partial^3 u}{\partial x^3} = 0$

### Boundary Condition Syntax
Boundary conditions must be passed as a semicolon-separated string of assignments:
```julia
"u(x=0, t) = 1.0; u(x=Inf, t) = 0"
```
Each condition defines the function, input restrictions, and its fixed value.

### Variable Decomposition
Variables are automatically split into:
- **Inputs** (e.g. `x`, `y`) — variables others depend on
- **Outputs** (e.g. `u(x,y)`) — variables defined as functions of inputs

The code then attempts similarity substitutions of the form:
$\eta = y x^m$ and  $u(x, y) = y^n f(\eta)$
and tries to find powers \( n, m \) that reduce the PDE.

---

## 🔬 Dilation Method (Recommended)

The dilation method finds similarity solutions by exact linear algebra — no guessing
or scanning. It assigns a formal scaling exponent to each variable, builds the
invariance system whose null space encodes all valid similarity transformations, and
derives the reduced ODE from the null-space vector.

### `find_ode_dilation` — Symbolic API

```julia
using SimilaritySolver, Symbolics

# Heat equation: u_t = u_xx
@variables x t u(..)
Dt = Differential(t); Dx = Differential(x)
heat = Dt(u(x,t)) - Dx(Dx(u(x,t)))
results = find_ode_dilation(heat; indep_vars=[x,t], dep_vars=[u(x,t)])

println(results[1]["similarity_variable"])  # x*t^(-1//2)
println(results[1]["gamma"])                # 0//1  (u scales as t^0)
println(results[1]["PDE"])                  # ODE: -(1//2)*η_bare*Differential(η_bare)(f_dil(...)) - Differential(η_bare)(Differential(η_bare)(f_dil(...)))
```

```julia
# KdV: u_t + 6u u_x + u_xxx = 0
@variables x t u(..)
Dt = Differential(t); Dx = Differential(x)
kdv = Dt(u(x,t)) + 6*u(x,t)*Dx(u(x,t)) + Dx(Dx(Dx(u(x,t))))
results = find_ode_dilation(kdv; indep_vars=[x,t], dep_vars=[u(x,t)])

println(results[1]["similarity_variable"])  # x*t^(-1//3)
println(results[1]["gamma"])                # -2//1  (u ~ t^(-2/3))
```

### `find_similarity_v2` — String API

`find_similarity_v2` wraps `find_ode_dilation` with the same string-based interface
as `find_similarity`:

```julia
result = find_similarity_v2(
    "du/dt - d2u/d2x = 0",
    "u(x=Inf, t) = 0; u(x, t=0) = 1"
)
println(result["similarity_variable"])  # η = x * t^(-1//2)
println(result["PDE_similarity"])       # the reduced heat ODE
```

### Result Dict Structure

Each entry in the `Vector{Dict}` returned by `find_ode_dilation` has the following keys:

| Key | Type | Description |
|-----|------|-------------|
| `"PDE"` | `Num` | Reduced ODE expression (set to zero) |
| `"PDE_similarity"` | `Num` | Same ODE, simplified |
| `"similarity_variable"` | `Num` | The similarity variable η in terms of original vars |
| `"output_similarity"` | `Dict` | Maps each dep. var to its ansatz (e.g. `u => t^γ·f(η)`) |
| `"gamma"` | `Rational` | Scaling exponent γ of the dependent variable |
| `"alpha_vec"` | `Vector{Rational}` | Normalized scaling exponents of the independent variables |
| `"gamma_vals"` | `Vector{Rational}` | Scaling exponent per dependent variable (systems) |

---

## 📖 API Overview

### `find_ode_dilation(pde::Num; indep_vars, dep_vars, verbose=false)`
Exact dilation-symmetry method. Recommended for all new use.

### `find_similarity_v2(pde::String, boundary_conditions::String; parameters=[], verbose=false)`
String-API wrapper for `find_ode_dilation`.

### `find_similarity(pde::String, boundary_conditions::String)`
Original heuristic solver (power-scanning). Kept for backwards compatibility.

### `parse_boundary_condition(condition::String)`
Parses a boundary condition string into symbolic structure.

### `parse_pde(expr::String, input_vars, output_vars)`
Parses symbolic PDE from string notation like `d2x/dy`.

### `decomposeVars(vars)`
Splits symbolic variables into inputs and outputs.

---

## 🧪 Testing

Run the full test suite:

```bash
julia --project test/runtests.jl
```

Expected: all `find_ode_dilation` tests pass; the old `find_ode` tests are marked
`Broken` due to a Symbolics.jl API change (see `SYMBOLICS_BUGS.md`).

---

## 🤝 Contributing

Contributions are welcome! Please open an issue or pull request if you want to:
- Add support for more complex boundary conditions
- Extend power-guessing logic
- Add LaTeX export or equation rendering

---

## 📄 License

MIT License. See `LICENSE` file for details.

---

## 📬 Contact

Maintained by [@elvispy](https://github.com/elvispy). For academic collaborations or bug reports, feel free to open an issue or email directly.
