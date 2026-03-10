# Comparison to Existing Tools

`SimilaritySolver.jl` occupies a unique niche in the landscape of symbolic PDE tools.

## Feature Matrix

| Feature | SimilaritySolver.jl | Maple (GeM/PDEtools) | Mathematica (MathLie) | SymPy (Python) |
|---|---|---|---|---|
| **Language** | Julia | Maple (Commercial) | Mathematica (Commercial) | Python (Open Source) |
| **Dilation Symmetry** | **Automated** | Manual/Semi-auto | Manual | Manual |
| **Coupled Systems** | **Yes** | Yes | Yes | Limited |
| **Exact Rational Exponents**| **Yes** | Yes | Yes | No (Heuristic) |
| **JOSS-Compliant Docs** | **Yes** | N/A | N/A | Yes |

## Research Impact and Statement of Need

As of 2026, similarity solutions remain a foundational technique in fluid mechanics, nonlinear diffusion, and materials science. However, these solutions are frequently still computed **by hand** or via commercial software like Maple or Mathematica. 

Existing open-source ecosystems (Julia, Python) lack a dedicated tool that takes an arbitrary PDE and automatically detects scaling (dilation) symmetries without manual input of the symmetry group.

`SimilaritySolver.jl` addresses this gap by providing:
1.  **Automation**: The invariance system is derived and solved via exact linear algebra—no "guessing" exponents.
2.  **Scalability**: Handles coupled systems of multiple PDEs, which are common in complex physical models.
3.  **Reproducibility**: Rational exponents ($1/2$, $1/3$) are preserved exactly, ensuring that downstream numerical solvers or PINNs receive citable, exact mathematical structures.

By providing a native Julia implementation, `SimilaritySolver.jl` enables researchers to integrate automated similarity reduction directly into modern scientific machine learning (SciML) workflows.
