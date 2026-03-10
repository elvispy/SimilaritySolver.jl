# Algorithm: Dilation Symmetry Method

`SimilaritySolver.jl` implements the **dilation (scaling) symmetry** method for PDE reduction. This approach is more systematic than traditional heuristic methods because it replaces trial-and-error with exact linear algebra.

## Core Concepts

The method assumes the PDE is invariant under a scaling transformation of the form:
$x_i \to \lambda^{a_i} x_i, \quad u_j \to \lambda^{c_j} u_j$

where $\lambda$ is a positive real scaling parameter. If the PDE is invariant under this transformation, we can find a **similarity variable** $\eta$ that remains unchanged when $\lambda$ scales.

## The 4-Step Process

### 1. Build the Invariance System

For each term in the PDE, we calculate its **scaling degree** based on the exponents $a_i$ and $c_j$. If the PDE consists of $N$ terms, all $N$ terms must scale by the same total degree for the equation to be invariant. This leads to a system of $N-1$ linear equations:
$\sum_{i} a_i \cdot (\text{indep\_deg}_{i,k}) + \sum_{j} c_j \cdot (\text{dep\_deg}_{j,k}) = \text{Constant}_k$

### 2. Solve the Null Space

We construct a matrix $A$ representing this system. The valid scaling exponents $(a_i, c_j)$ are the vectors in the **null space** of $A$. `SimilaritySolver.jl` uses **exact rational arithmetic** to find this null space, ensuring that exponents like $1/2$ or $2/3$ are identified exactly rather than as floating-point approximations.

### 3. Generate Similarity Candidates

Each null-space vector defines a potential similarity transformation. We construct the similarity variable $\eta$ and the dependent-variable ansatz:
- **Similarity Variable**: $\eta = x_1 \cdot x_2^{-a_1/a_2}$ (choosing a "pivot" variable $x_2$).
- **Ansatz**: $u = x_2^{c/a_2} \cdot f(\eta)$.

### 4. PDE to ODE Reduction

The final step substitutes the ansatz into the original PDE. Using the chain rule, all occurrences of the original independent variables $(x_i, t)$ must cancel out, leaving a pure ODE in terms of $f(\eta)$ and its derivatives. 

`SimilaritySolver.jl` performs this substitution symbolically and simplifies the result to produce a clean, machine-readable ODE.
