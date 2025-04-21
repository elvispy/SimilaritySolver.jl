[![Build Status](https://github.com/elvispy/SimilaritySolver.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/elvispy/SimilaritySolver.jl/actions/workflows/CI.yml?query=branch%3Amain)

# Symbolic Similarity PDE Solver (Julia)

This Julia package provides a set of tools for symbolic manipulation and similarity transformation of partial differential equations (PDEs). It is built on top of the `Symbolics.jl` and `SymbolicUtils.jl` ecosystems and supports automatic classification of variables, similarity analysis, boundary condition parsing, and ODE reduction.

---

## 🚀 Features

- **Automatic decomposition** of symbolic variables into inputs and outputs.
- **Similarity transformation detection** for 2D PDEs.
- **Boundary condition parsing** and transformation under similarity variables.
- **Power guessing heuristic** for identifying simplifying substitutions.
- **Parsing PDEs from strings** using symbolic differential notation like `d2x/dy`.

---

## 📦 Installation

Make sure you have Julia installed. Then clone this repository and instantiate dependencies:

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

## 📖 API Overview

### `find_ode(symbolicPDE::Num; vars::Vector{Num})`
Attempts to find a similarity transformation and reduce the given PDE to an ODE.

### `find_similarity(pde::String, boundary_conditions::String)`
Full pipeline: parse PDE + boundary conditions, attempt similarity reduction.

### `parse_boundary_condition(condition::String)`
Parses a boundary condition string into symbolic structure.

### `parse_pde(expr::String, input_vars, output_vars)`
Parses symbolic PDE from string notation like `d2x/dy`.

### `decomposeVars(vars)`
Splits symbolic variables into inputs and outputs.

---

## 🧪 Testing

You can run the test examples at the bottom of the main file to validate the behavior. For full test coverage, consider structuring tests using `Test.jl` in a `test/` folder.

---

## 🧠 Motivation

Similarity solutions help reduce PDEs to ODEs, making complex problems more tractable. This symbolic toolchain is meant to facilitate experimentation with such reductions, automating variable classification and heuristic power substitution.

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

