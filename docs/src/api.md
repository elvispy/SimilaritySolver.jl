# API Reference

`SimilaritySolver.jl` exports three primary functions for similarity reduction.

## Symbolic API

`find_ode_dilation` is the core solver. It uses exact scaling-symmetry (dilation) analysis and is recommended for most applications.

```@docs
find_ode_dilation
```

## String API (High-Level)

`find_similarity_v2` provides a high-level string-based interface that parses PDEs and boundary conditions automatically.

```@docs
find_similarity_v2
```

## Heuristic Solver (Legacy)

`find_similarity` is the original heuristic solver that scans power-law combinations. It is kept for backward compatibility.

```@docs
find_similarity
```
