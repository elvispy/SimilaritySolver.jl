---
title: 'SimilaritySolver.jl: Automated Scaling Symmetry Detection and PDE Reduction in Julia'
tags:
  - Julia
  - Symbolics
  - Partial Differential Equations
  - Similarity Solutions
  - Dilation Symmetry
authors:
  - name: Elvis Aguero
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 10 March 2026
bibliography: paper.bib
---

# Summary

Finding similarity solutions for Partial Differential Equations (PDEs) is a fundamental technique in mathematical physics and fluid mechanics. Traditionally, this process requires manually identifying scaling invariances (dilation symmetries) by guessing exponents that balance terms within the equation. `SimilaritySolver.jl` automates this entire pipeline. It takes an arbitrary symbolic PDE or coupled system, detects valid dilation symmetries using exact rational linear algebra, and derives the reduced Ordinary Differential Equation (ODE) system automatically.

# Statement of Need

As of 2026, many similarity solutions in fluid mechanics, nonlinear diffusion, and materials science are still computed manually or via proprietary commercial software like Maple's GeM package [@cheviakov2007gem] or Mathematica's MathLie. Existing open-source symbolic ecosystems, such as SymPy in Python or the core `Symbolics.jl` in Julia, lack a dedicated, high-level pipeline for automated scaling-symmetry detection.

The need for automated symmetry analysis has grown with the rise of scientific machine learning (SciML). Researchers training Physics-Informed Neural Networks (PINNs) often compute Lie symmetries manually to validate their models or regularize training [@eshima2025similarity]. `SimilaritySolver.jl` addresses this gap by providing a machine-readable API that integrates directly with the Julia SciML ecosystem [@gowda2021high]. By automating the derivation of exact similarity variables, the package enables researchers to focus on the physical interpretation of invariant solutions rather than the algebraic drudgery of their derivation.

# Algorithm and Implementation

`SimilaritySolver.jl` implements the dilation symmetry method. The software constructs an invariance system by calculating the scaling degree of each term in the PDE. This results in a system of linear equations where the valid scaling exponents reside in the null space of the invariance matrix.

A key architectural feature of the package is its use of exact rational arithmetic (`Rational{Int}`). This ensures that scaling exponents like $1/2$ or $2/3$ are identified exactly, avoiding the numerical drift and heuristic "snapping" common in floating-point implementations. Furthermore, the package includes a principled workaround for a known limitation in the underlying symbolic engine's handling of high-order Leibniz expansions, ensuring mathematical correctness for PDEs of any order.

# Functionality and Examples

The package supports both a high-level string-based API and a full symbolic API built on `Symbolics.jl`.

## Heat Equation
For the heat equation $u_t = u_{xx}$, the solver automatically identifies the classic similarity variable $\eta = x/\sqrt{t}$ and the reduced ODE $f'' + \frac{1}{2}\eta f' = 0$.

## Coupled Systems
`SimilaritySolver.jl` can handle complex, multi-variable systems. We validated the package against the Inertio-Marangoni (IM) regime of inertial surfactant dynamics from @eshima2025similarity. The solver correctly recovered the shared dilation symmetry for the coupled mass, momentum, and surfactant transport equations, producing three coupled ODEs free of the original independent variables.

# Comparison to Existing Work

Unlike the recent "Invariant Reduction" frameworks implemented in Maple [@druzhkov2024invariant; @druzhkov2025invariant], which require users to provide the symmetry group manually, `SimilaritySolver.jl` automatically detects dilation symmetries without user intervention. While it is more specialized than general Lie-symmetry packages, its focus on scaling symmetries allows for a highly automated and robust workflow tailored to the most common symmetry found in physical laws.

# Acknowledgements

The author would like to thank the Julia Symbolics community for their foundational work.

# References
