# Symbolics.jl / SymbolicUtils.jl Bug Tracker

Bugs discovered while developing `SimilaritySolver.jl`.  Each entry has a
minimal reproduction script and the observed vs. expected behaviour.

---

## Bug 1 — Incorrect Leibniz coefficients in `expand_derivatives` for high-order derivatives of `x^γ · f(η(x,t))`

**Packages**: `Symbolics.jl` (confirmed on v6.22.0), `SymbolicUtils.jl` (v3.7.2)
**Discovered**: 2026-03 during KdV similarity reduction in SimilaritySolver.jl

### Summary

`expand_derivatives(Differential(x)^3(x^(-2//1) * f_fn(η_sym)))` returns a
coefficient of **-4** for the `f''(η)·(∂η/∂x)²` term.  The correct Leibniz
formula gives **-6** (= 3·C(3,1)·D_x(x⁻²) = 3·(-2x⁻³)).

This error also manifests for any rational (or negative integer) exponent γ and
any order n ≥ 3 where the mixed Leibniz term
`C(n,k)·D_x^k(x^γ)·D_x^{n-k}(f(η))` is affected.

### Minimal reproduction

```julia
using Symbolics, SymbolicUtils
@variables x t

η_fn = Symbolics.variable(:η, T=Symbolics.FnType{Tuple{Any,Any}, Real})
η_sym = η_fn(x, t)
f_fn  = Symbolics.variable(:f, T=Symbolics.FnType{Tuple{Any}, Real})

# Build chain subs for KdV: η = x * t^(-1/3)
η_expr = x * t^(-1//3)
chain_subs = Dict{Any,Any}(
    Differential(x)(η_sym) => t^(-1//3),
    Differential(t)(η_sym) => -x/(3*t^(4//3)),
    Differential(x)(Differential(x)(η_sym)) => 0,
)

# === BUGGY path ===
e_bug = expand_derivatives(Differential(x)(Differential(x)(Differential(x)(x^(-2//1) * f_fn(η_sym)))))
e_bug2 = substitute(e_bug, chain_subs)
# ... simplify and substitute t → x^3/η^3, x → 1
# Yields coefficient -4 for f''(η)·η² term

# === CORRECT path (manual Leibniz) ===
# D_x^3(x^{-2} · f(η)) = Σ_{k=0}^3 C(3,k) · [−2]↓k · x^{-2-k} · D_x^{3-k}(f(η))
# k=1 term: C(3,1)·(-2)·x^{-3}·D_x^2(f(η)) = -6x^{-3}·f''·(D_x η)²
# This gives coefficient -6 as expected.
```

### Root cause (hypothesis)

When Symbolics internally represents `x^{-2}` via `exp(-2·log(x))` for the
purpose of differentiation, the repeated product-rule application in
`expand_derivatives` leaves `Differential(x)(x)` unevaluated in some
intermediate terms.  These residual `Differential(x)(x)` tokens are later
evaluated to 1, but not before certain `D_x^k(x^γ)·D_x^{n-k}(f)` cross-terms
have already been counted twice (while one copy that should be included from
`D_x(x^γ)·D_x^{n-1}(f)` in the last Leibniz step is lost).

### Workaround (used in SimilaritySolver.jl)

Apply the Leibniz rule manually:

```julia
function _ffact(γ::Rational{Int}, k::Int)::Rational{Int}
    k == 0 && return 1//1
    prod(γ - i for i in 0:k-1)
end

function manual_Dx_n(pivot, γ, f_fn, η_sym, n)
    result = Symbolics.Num(0)
    for k in 0:n
        ff = _ffact(γ, k)
        iszero(ff) && continue
        u_k = ff * pivot^(γ - k)
        v = f_fn(η_sym)
        for _ in 1:(n-k); v = Differential(pivot)(v); end
        v = expand_derivatives(v)   # only the PURE f-composition, no x^γ
        result += binomial(n, k) * u_k * v
    end
    return result
end
```

`expand_derivatives` applied to the *pure-composition* part `D_x^j(f(η(x,t)))`
works correctly; only the joint product `D_x^n(x^γ · f(η))` is affected.

### Impact

Any PDE reduction that involves substituting a power-law ansatz
`u = x^γ · f(η)` into a ≥ 3rd-order derivative of `u` will produce **wrong
ODE coefficients** if `expand_derivatives` is used directly on the substituted
expression.

---
