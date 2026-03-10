# Architectural Decisions and Trade-offs

This document outlines the key design decisions made during the development of `SimilaritySolver.jl` and the trade-offs involved in its implementation.

## 1. Exact Rational Arithmetic vs. Floating-Point

The most significant architectural choice in `SimilaritySolver.jl` is the use of **exact rational arithmetic** for all scaling-symmetry computations.

### Decision
The invariance system and its null space are solved using `Rational{Int}` rather than `Float64`.

### Rationale
- **Mathematical Rigor**: Similarity exponents in physics are almost always simple rationals ($1/2$, $1/3$, $2/5$). Numerical solvers often return values like $0.499999998$, which require heuristic "snapping" to the nearest fraction.
- **Robustness**: Exact arithmetic prevents numerical drift during the substitution phase, where high-order derivatives of terms like $x^{1/2}$ must cancel out precisely.

### Trade-off
- **Performance**: Rational arithmetic is computationally more expensive than floating-point. However, since the invariance systems for most PDEs are small (typically $< 10 \times 10$ matrices), the absolute overhead is negligible compared to the symbolic manipulation time.

## 2. Symbolic Leibniz Workaround

During development, we identified a limitation in the underlying symbolic engine (`Symbolics.jl`) regarding the expansion of high-order derivatives of product terms.

### Decision
We implemented a manual, principled application of the Leibniz rule for the dependent-variable ansatz.

### Rationale
Directly calling `expand_derivatives` on $D_x^n(x^\gamma \cdot f(\eta))$ can produce incorrect coefficients for orders $n \geq 3$ due to internal simplification orderings. By manually applying the Leibniz rule, we ensure the mathematical correctness of the reduced ODE regardless of upstream bugs.

## 3. Dependent Variable Extraction Strategy

### Decision
In the high-level String API (`find_similarity_v2`), dependent variables are extracted directly from the PDE string using regex, rather than being inferred from the boundary conditions.

### Rationale
- **Parsing Robustness**: Boundary conditions often involve derivatives (e.g., $d\psi/dy(x, y=0) = 0$). Inferring variables from these strings can lead to misidentifying "dy" as a new dependent variable, which poisons the invariance matrix and crashes the solver.
- **Single Source of Truth**: The PDE itself is the definitive source for which variables are dependent.

## 4. Parsimony-Based Null Space Expansion

### Decision
When the dilation null space is multi-dimensional (e.g., in the Blasius equation), the solver enumerates integer linear combinations of basis vectors and sorts them by their $L_1$ norm.

### Rationale
The "physically interesting" similarity variable is often a simple combination (like $v_1 + v_2$) rather than a raw basis vector. Sorting by the $L_1$ norm (Parsimony Ordering) ensures that the most likely physical candidates are presented to the user first.

---

*These decisions were made to prioritize the mathematical reliability and user-friendliness of SimilaritySolver.jl for the research community.*
