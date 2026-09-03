# Known Limitations

While `SimilaritySolver.jl` is a powerful tool for PDE reduction, users should be aware of the following limitations.

## Nonlocal operators

Fourier multipliers of fixed homogeneity are supported through the `nonlocal_ops` keyword:

```julia
@variables x t q(..) u(..) H(..)
# Hilbert transform: degree 0 (scale invariant)
find_dilation_symmetry(pdes, [x,t], [q(x,t), u(x,t)]; nonlocal_ops = Dict(:H => 0))
# fractional Laplacian |k|^α dual to x
find_dilation_symmetry(eq,  [x,t], [u(x,t)];          nonlocal_ops = Dict(:Λ => (:x => 3//2)))
```

This affects **symmetry detection only**. A nonlocal PDE reduces to an integro-differential
equation rather than an ODE, so `reduce_to_ode` does not handle these systems; use
`find_dilation_symmetry` / `build_invariance_system` and carry out the reduction by hand.

## Scope of Symmetries

`SimilaritySolver.jl` is specifically designed to find **dilation (scaling) symmetries**. It does not currently support:
- **Translational Symmetries**: $x \to x + c$ (Traveling wave solutions).
- **Rotational Symmetries**.
- **General Lie Symmetries**: Full infinitesimal generator analysis.

For these more general cases, classical Lie-symmetry software (often in Maple or Mathematica) is required.

## Symbolics.jl Bug (Workaround Applied)

We discovered a bug in `Symbolics.jl` / `SymbolicUtils.jl` where `expand_derivatives` produces incorrect Leibniz coefficients for high-order derivatives ($n \geq 3$) of a product like $x^\gamma \cdot f(\eta)$. 

**Example**: $D_x^3(x^{-2} \cdot f(x))$ incorrectly returns a coefficient of $-4$ for the $f''(x)$ term instead of the correct $-6$.

`SimilaritySolver.jl` includes a **principled workaround** that manually applies the Leibniz rule for the dependent-variable ansatz, ensuring that all reduced ODEs are mathematically correct even when the underlying symbolic engine has this bug.

## Dimensionality

- **Null Space Complexity**: For PDEs with very high-dimensional null spaces, the number of potential similarity candidates can grow. The solver uses **parsimony ordering** to return the simplest combinations first.
- **BC Transformation**: Full boundary condition transformation for the new dilation method is currently under development (BCs are parsed but not yet fully transformed into the similarity domain in `find_similarity_v2`).
