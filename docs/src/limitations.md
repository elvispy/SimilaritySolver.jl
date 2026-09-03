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

## Exact vs. asymptotic invariance

`find_dilation_symmetry` tests **exact** dilation invariance. Most published models are not
scale invariant as written — a closure (a logarithmic isotherm, an entropic free energy, a grain
diameter, a Hamaker constant) pins a length. Physical similarity solutions usually live in an
*asymptotic* regime where those terms are subdominant, which is Barenblatt's intermediate
asymptotics; a scale may be present and a similarity solution still exist.

Use [`find_dominant_balances`](@ref) for that case. It enumerates sub-balances, keeps those with
a dilation symmetry, and — crucially — applies the self-consistency test that every dropped term
really is negligible under the scaling the retained terms imply. On the depth-averaged granular
model of Deléage & Richard (JFM 1009 A57, 2025) the full system has nullity 0, raw enumeration
returns 17,643 homogeneous subsets, and the self-consistency filter leaves 2.

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
