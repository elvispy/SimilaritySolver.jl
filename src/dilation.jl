using LinearAlgebra

# Module-level similarity variable (used as placeholder across all reductions)
@variables η_bare f_bare

# ---------------------------------------------------------------------------
# Scaling-degree computation
# ---------------------------------------------------------------------------
# Each variable is assigned a formal exponent parameter (a Symbol).
# The scaling degree of an expression is a Dict{Symbol,Rational{Int}}
# representing the linear combination of those parameter Symbols.

const ScalingDegree = Dict{Symbol, Rational{Int}}

_add_degrees(d1::ScalingDegree, d2::ScalingDegree) = mergewith(+, d1, d2)
_scale_degree(d::ScalingDegree, s::Rational{Int})  = Dict(k => s * v for (k, v) in d)
_neg_degree(d::ScalingDegree)                       = _scale_degree(d, -1//1)
_zero_degree()                                       = ScalingDegree()

"""
    _make_maps(indep_vars, dep_vars)
    -> (indep_map::Dict{Symbol,Symbol}, dep_map::Dict{Symbol,Symbol})

Build name-based lookup tables used by `_scaling_degree_raw`:
- `indep_map`: :x => :a_x  for each independent variable x
- `dep_map`:   :u => :c_u  for each dependent variable u (stored as u(x,...))
"""
function _make_maps(indep_vars, dep_vars)
    indep_map = Dict{Symbol,Symbol}()
    for v in indep_vars
        raw = Symbolics.unwrap(v)
        nm  = SymbolicUtils.issym(raw) ? raw.name : Symbol(v)
        indep_map[nm] = Symbol("a_", nm)
    end
    dep_map = Dict{Symbol,Symbol}()
    for v in dep_vars
        raw = Symbolics.unwrap(v)
        # u(x,t) is a Term whose operation is Sym{FnType}
        nm = if SymbolicUtils.iscall(raw)
            SymbolicUtils.operation(raw).name
        else
            SymbolicUtils.issym(raw) ? raw.name : Symbol(v)
        end
        dep_map[nm] = Symbol("c_", nm)
    end
    return indep_map, dep_map
end

function _param_names(indep_vars, dep_vars)
    indep_map, dep_map = _make_maps(indep_vars, dep_vars)
    return vcat(collect(values(indep_map)), collect(values(dep_map))),
           indep_map, dep_map
end

"""
    _scaling_degree_raw(raw, indep_map, dep_map) -> ScalingDegree

Recursively compute the ε-scaling degree of a raw SymbolicUtils expression.
"""
function _scaling_degree_raw(raw, indep_map::Dict{Symbol,Symbol}, dep_map::Dict{Symbol,Symbol})
    # --- concrete numbers ---
    raw isa Number && return _zero_degree()

    # --- atomic Sym: independent variable or parameter ---
    if SymbolicUtils.issym(raw)
        nm = raw.name
        if haskey(indep_map, nm)
            return ScalingDegree(indep_map[nm] => 1//1)
        end
        return _zero_degree()   # parameter
    end

    SymbolicUtils.iscall(raw) || return _zero_degree()

    op   = SymbolicUtils.operation(raw)
    args = SymbolicUtils.arguments(raw)

    # --- Differential(xᵢ)(inner) ---
    if op isa Differential
        inner_deg = _scaling_degree_raw(args[1], indep_map, dep_map)
        dx_raw = Symbolics.unwrap(op.x)
        dx_nm  = SymbolicUtils.issym(dx_raw) ? dx_raw.name : nothing
        if dx_nm !== nothing && haskey(indep_map, dx_nm)
            return _add_degrees(inner_deg, ScalingDegree(indep_map[dx_nm] => -1//1))
        end
        return inner_deg
    end

    # --- Dependent variable call: u(x, t) ---
    if SymbolicUtils.issym(op)
        nm = op.name
        if haskey(dep_map, nm)
            return ScalingDegree(dep_map[nm] => 1//1)
        end
        return _zero_degree()
    end

    # --- Addition ---
    if op === (+)
        degrees = [_scaling_degree_raw(a, indep_map, dep_map) for a in args]
        non_zero = filter(d -> !isempty(d), degrees)
        isempty(non_zero) && return _zero_degree()
        return first(non_zero)
    end

    # --- Multiplication ---
    if op === (*)
        return reduce(_add_degrees,
                      [_scaling_degree_raw(a, indep_map, dep_map) for a in args];
                      init=_zero_degree())
    end

    # --- Power: Pow(base, exp) ---
    if op === (^)
        base_raw, exp_raw = args
        exp_val = try
            Float64(Symbolics.unwrap(Symbolics.Num(exp_raw)))
        catch
            return _zero_degree()
        end
        isfinite(exp_val) || return _zero_degree()
        exp_rat = rationalize(Int, exp_val; tol=1e-9)
        return _scale_degree(_scaling_degree_raw(base_raw, indep_map, dep_map), exp_rat)
    end

    # --- Division ---
    if op === (/)
        num_deg = _scaling_degree_raw(args[1], indep_map, dep_map)
        den_deg = _scaling_degree_raw(args[2], indep_map, dep_map)
        return _add_degrees(num_deg, _neg_degree(den_deg))
    end

    # Fallback
    return _zero_degree()
end


# ---------------------------------------------------------------------------
# Build the invariance linear system
# ---------------------------------------------------------------------------

"""
    _collect_add_terms(raw) -> Vector

Return the top-level additive terms of `raw`.
"""
function _collect_add_terms(raw)
    if SymbolicUtils.iscall(raw) && SymbolicUtils.operation(raw) === (+)
        return SymbolicUtils.arguments(raw)
    end
    return [raw]
end

"""
    build_invariance_system(pde, indep_vars, dep_vars)
    -> (A::Matrix{Float64}, param_names::Vector{Symbol})

Matrix A of size (n_terms-1, n_params). Its null space gives the valid
dilation exponents. Columns are ordered as [a_x1, …, a_xn, c_u1, …, c_um].
"""
function build_invariance_system(pde, indep_vars::Vector, dep_vars::Vector)
    indep_map, dep_map = _make_maps(indep_vars, dep_vars)
    # stable column ordering
    param_names = vcat([indep_map[k] for k in [Symbolics.unwrap(v).name for v in indep_vars]],
                       [dep_map[begin_dep_name(v)] for v in dep_vars])

    raw   = Symbolics.unwrap(pde)
    terms = _collect_add_terms(raw)
    degs  = [_scaling_degree_raw(t, indep_map, dep_map) for t in terms]

    n_params = length(param_names)
    n_terms  = length(degs)

    n_terms <= 1 && return zeros(Float64, 0, n_params), param_names

    A = zeros(Float64, n_terms - 1, n_params)
    ref = degs[1]
    for i in 2:n_terms
        diff = _add_degrees(degs[i], _neg_degree(ref))
        for (j, pname) in enumerate(param_names)
            A[i-1, j] = Float64(get(diff, pname, 0//1))
        end
    end
    return A, param_names
end

# Helper: get the canonical name used in dep_map
function begin_dep_name(v)
    raw = Symbolics.unwrap(v)
    if SymbolicUtils.iscall(raw)
        return SymbolicUtils.operation(raw).name
    end
    return SymbolicUtils.issym(raw) ? raw.name : Symbol(v)
end

# Helper: get canonical name for indep var
function begin_indep_name(v)
    raw = Symbolics.unwrap(v)
    return SymbolicUtils.issym(raw) ? raw.name : Symbol(v)
end


# ---------------------------------------------------------------------------
# Find dilation symmetry (null space)
# ---------------------------------------------------------------------------

"""
    find_dilation_symmetry(pde, indep_vars, dep_vars)
    -> (cols::Vector{Vector{Rational{Int}}}, param_names::Vector{Symbol})

Return null-space vectors of the invariance system as rational vectors,
computed exactly via Gaussian elimination over ℚ.
"""
function find_dilation_symmetry(pde, indep_vars::Vector, dep_vars::Vector)
    A, param_names = build_invariance_system(pde, indep_vars, dep_vars)
    n = length(param_names)

    if size(A, 1) == 0
        return [Vector{Rational{Int}}([i == k ? 1//1 : 0//1 for i in 1:n]) for k in 1:n],
               param_names
    end

    # Exact rational null space via row reduction
    cols = _rational_nullspace(A, n)
    return cols, param_names
end

"""
    _rational_nullspace(A::Matrix{Float64}, n::Int) -> Vector{Vector{Rational{Int}}}

Compute the null space of A exactly by converting to rational, performing
row reduction over ℚ, and reading off the free variables.
"""
function _rational_nullspace(A::Matrix{Float64}, n::Int)
    m = size(A, 1)
    # Round to nearest integer matrix (the invariance constraints are always integers)
    Aint = round.(Int, A)
    # Work over Rational{Int64} to stay exact
    Arat = Rational{Int64}.(Aint)

    # Augment with identity to track column combinations: [A | I]
    aug = hcat(Arat, Matrix{Rational{Int64}}(I, m, m))

    # Row-echelon form over ℚ (partial pivoting)
    pivot_cols = Int[]
    row = 1
    for col in 1:n
        # Find pivot
        piv = findfirst(r -> !iszero(aug[r, col]), row:m)
        piv === nothing && continue
        piv += row - 1
        # Swap
        aug[row, :], aug[piv, :] = aug[piv, :], aug[row, :]
        # Eliminate
        aug[row, :] ./= aug[row, col]
        for r in 1:m
            r == row && continue
            aug[r, :] .-= aug[r, col] .* aug[row, :]
        end
        push!(pivot_cols, col)
        row += 1
    end

    free_cols = setdiff(1:n, pivot_cols)
    isempty(free_cols) && return Vector{Rational{Int}}[]

    # Build null space vectors: for each free variable j, set it to 1, others to 0,
    # solve for pivot variables.
    null_vecs = Vector{Rational{Int}}[]
    for j in free_cols
        v = zeros(Rational{Int64}, n)
        v[j] = 1//1
        for (k, pc) in enumerate(pivot_cols)
            # aug[k, pc] == 1 after RREF; aug[k, j] gives the coefficient
            v[pc] = -aug[k, j]
        end
        # Normalize: make first nonzero entry positive, clear denominators
        first_nz = findfirst(!iszero, v)
        if first_nz !== nothing && v[first_nz] < 0
            v = -v
        end
        # Scale to integer by clearing denominators
        denom = lcm([x.den for x in v if !iszero(x)])
        v_int = [Rational{Int64}(numerator(x * denom), 1) for x in v]
        push!(null_vecs, v_int)
    end
    return null_vecs
end


# ---------------------------------------------------------------------------
# Build similarity ansatz from a null-space vector
# ---------------------------------------------------------------------------

"""
    build_similarity_candidates(null_vec, param_names, indep_vars, dep_vars)
    -> Vector of NamedTuples

For each non-zero indep exponent as pivot, produce one candidate:
  η = pivot * ∏ other_xᵢ^{-αᵢ/α_pivot}
  u = pivot^{γ/α_pivot} * f(η)
"""
function build_similarity_candidates(null_vec, param_names, indep_vars, dep_vars)
    n_i = length(indep_vars)
    α   = null_vec[1:n_i]
    γ   = null_vec[n_i+1:end]

    results = Vector{Any}()
    for (pivot_idx, αp) in enumerate(α)
        iszero(αp) && continue
        α_norm = α .// αp
        γ_norm = γ .// αp

        xp = indep_vars[pivot_idx]
        η_expr = xp
        for i in 1:n_i
            i == pivot_idx && continue
            iszero(α_norm[i]) && continue
            # invariant combination: x_pivot * x_i^{-a_pivot/a_i}
            # with α_norm[i] = a_i/a_pivot, the exponent is -1/α_norm[i]
            exp_i = -(α_norm[i].den // α_norm[i].num)
            η_expr = η_expr * indep_vars[i]^exp_i
        end

        for (j, dep) in enumerate(dep_vars)
            push!(results, (
                η_expr     = η_expr,
                pivot      = xp,
                pivot_idx  = pivot_idx,
                gamma      = γ_norm[j],
                dep_idx    = j,
                alpha_norm = α_norm,
                gamma_norm = γ_norm,
            ))
        end
    end
    return results
end


# ---------------------------------------------------------------------------
# Reduce PDE to ODE
# ---------------------------------------------------------------------------

# ── Leibniz helpers (bypass Symbolics bug in D_x^n(x^γ·f)) ─────────────────
# Symbolics.expand_derivatives incorrectly computes D_x^n(x^γ · f(η)) for
# n ≥ 3 and rational γ ≤ -1: the Leibniz coefficient for certain f'' terms
# is off by one.  We work around this by applying the Leibniz product rule
# *manually* and calling expand_derivatives only on each pure-f factor.

"""
    _ffact(γ, k) -> Rational{Int}
Falling factorial γ·(γ-1)·…·(γ-k+1).
"""
function _ffact(γ::Rational{Int}, k::Int)::Rational{Int}
    k == 0 && return 1//1
    prod(γ - i for i in 0:k-1)
end

"""
    _leibniz_dep_val(indep_vars, pivot_idx, gamma, f_fn, η_sym, orders)

Compute  ∂^{orders[1]}_{x1} · ∂^{orders[2]}_{x2} · … · (pivot^γ · f(η_sym))
using the Leibniz rule *only* on the pivot variable (for which pivot^γ has a
non-trivial derivative).  All other variables differentiate only through f(η_sym)
because ∂_xi(pivot) = 0 for xi ≠ pivot.
"""
function _leibniz_dep_val(indep_vars::Vector, pivot_idx::Int,
                           gamma::Rational{Int}, f_fn, η_sym,
                           orders::Vector{Int})
    pivot_var = indep_vars[pivot_idx]
    k_piv     = orders[pivot_idx]
    result    = Symbolics.Num(0)
    for j in 0:k_piv
        ff = _ffact(gamma, j)
        iszero(ff) && continue
        u_j = ff * pivot_var^(gamma - j)
        # Build D_pivot^(k_piv-j) · Π_{i≠pivot} D_xi^orders[i] of f(η_sym)
        v = f_fn(η_sym)
        for i in eachindex(indep_vars)
            cnt = (i == pivot_idx) ? (k_piv - j) : orders[i]
            for _ in 1:cnt; v = Differential(indep_vars[i])(v); end
        end
        v      = expand_derivatives(v)
        result = result + binomial(k_piv, j) * u_j * v
    end
    return result
end

"""
    _scan_dep_keys(raw, dep_raw) -> Vector{Any}

Walk the expression tree rooted at `raw` and collect every sub-expression
that is a derivative chain applied directly to `dep_raw` (including `dep_raw`
itself).  Returns the raw (unwrapped) expressions for use as substitution keys.
"""
function _scan_dep_keys(raw, dep_raw)
    keys = Any[]
    # Base: is this node dep_raw itself?
    if isequal(Symbolics.Num(raw), Symbolics.Num(dep_raw))
        push!(keys, raw)
        return keys
    end
    SymbolicUtils.iscall(raw) || return keys
    op   = SymbolicUtils.operation(raw)
    args = SymbolicUtils.arguments(raw)
    if op isa Differential
        # This node is D_xi(inner); check if inner is a dep-chain
        inner_keys = _scan_dep_keys(args[1], dep_raw)
        if !isempty(inner_keys)
            push!(keys, raw)   # this outer node is also a dep-chain key
        end
    end
    # Always recurse into all children (dep may appear nested inside products, etc.)
    for a in args
        append!(keys, _scan_dep_keys(a, dep_raw))
    end
    return keys
end

"""
    _dep_diff_orders(key_raw, dep_raw, indep_vars) -> Vector{Int}

Given a dep-chain key (the raw form of D_xi^k1(D_xj^k2(...(dep)...))),
return the differentiation order vector `orders` such that `orders[i]` is the
total number of times `indep_vars[i]` was differentiated.
"""
function _dep_diff_orders(key_raw, dep_raw, indep_vars)
    orders = zeros(Int, length(indep_vars))
    cur = key_raw
    while !isequal(Symbolics.Num(cur), Symbolics.Num(dep_raw))
        SymbolicUtils.iscall(cur) || break
        op = SymbolicUtils.operation(cur)
        op isa Differential || break
        xi_raw = Symbolics.unwrap(op.x)
        for (i, v) in enumerate(indep_vars)
            if isequal(Symbolics.Num(xi_raw), Symbolics.Num(Symbolics.unwrap(v)))
                orders[i] += 1
                break
            end
        end
        cur = SymbolicUtils.arguments(cur)[1]
    end
    return orders
end

"""
    _build_dep_subs(pde, dep, indep_vars, pivot_idx, gamma, f_fn, η_sym)

Scan `pde` for every derivative pattern applied to `dep`, then return a
substitution dictionary mapping each pattern to the correct Leibniz expansion
of  pivot^γ · f(η_sym).

This is order-agnostic: it works for PDEs of any derivative order without any
pre-specified maximum.
"""
function _build_dep_subs(pde, dep, indep_vars::Vector, pivot_idx::Int,
                          gamma::Rational{Int}, f_fn, η_sym)
    dep_raw = Symbolics.unwrap(dep)
    pde_raw = Symbolics.unwrap(pde)

    # Collect all dep-chain keys from the PDE tree (deduplicated)
    raw_keys = unique(_scan_dep_keys(pde_raw, dep_raw))

    subs = Dict{Any,Any}()
    for rk in raw_keys
        orders = _dep_diff_orders(rk, dep_raw, indep_vars)
        val    = _leibniz_dep_val(indep_vars, pivot_idx, gamma, f_fn, η_sym, orders)
        subs[Symbolics.Num(rk)] = val
    end
    return subs
end

"""
    reduce_to_ode(pde, indep_vars, dep_vars, η_expr, pivot, pivot_idx, gamma, dep_idx, alpha_norm)
    -> Num or nothing

Substitute the ansatz `dep = pivot^γ · f(η_sym)` into `pde`, apply chain-rule
substitutions, eliminate non-pivot variables via the similarity relation, then
divide out the pivot power and check no independent variables remain.

NOTE: step 1 uses a manual Leibniz expansion (not expand_derivatives on the full
product) to work around a Symbolics.jl bug where D_x^n(x^γ·f(η)) gives wrong
Leibniz coefficients for n ≥ 3, rational γ < 0.
"""
function reduce_to_ode(pde, indep_vars::Vector, dep_vars::Vector,
                        η_expr, pivot, pivot_idx::Int,
                        gamma::Rational{Int}, dep_idx::Int,
                        alpha_norm::Vector{Rational{Int}})
    dep = dep_vars[dep_idx]
    S(e, d) = simplify(substitute(e, d); expand=true)

    # Define η as a FUNCTION of the independent variables so that chain rule works.
    n_args    = length(indep_vars)
    any_tuple = Tuple{[Any for _ in 1:n_args]...}
    η_fn  = Symbolics.variable(:η_dil, T=Symbolics.FnType{any_tuple, Real})
    η_sym = η_fn(indep_vars...)
    f_fn  = Symbolics.variable(:f_dil, T=Symbolics.FnType{Tuple{Any}, Real})

    # 1. Scan pde for every derivative pattern on dep, build substitution dict
    #    via the Leibniz rule.  Order-agnostic: works for any PDE order.
    #    Avoids the Symbolics bug in expand_derivatives(D_x^n(x^γ·f(η))).
    dep_subs = _build_dep_subs(pde, dep, indep_vars, pivot_idx, gamma, f_fn, η_sym)
    expr = substitute(pde, dep_subs)

    # 2. Substitute chain-rule derivatives of η_sym (up to order 4) with concrete values.
    chain_subs = Dict{Any,Any}()
    function _add_chain_subs!(subs, sym_key, expr_val, depth)
        depth == 0 && return
        for xi in indep_vars
            new_key = Differential(xi)(sym_key)
            new_val = expand_derivatives(Differential(xi)(expr_val))
            subs[new_key] = new_val
            _add_chain_subs!(subs, new_key, new_val, depth - 1)
        end
    end
    _add_chain_subs!(chain_subs, η_sym, η_expr, 4)
    expr = S(expr, chain_subs)

    # 3. Hide η_sym → η_bare, then expand_derivatives to zero out D_xi(η_bare) = 0.
    expr = substitute(expr, η_sym => η_bare)
    expr = simplify(expand_derivatives(expr); expand=true)

    # 4. Express every non-pivot variable through η_bare:
    #    xi = pivot^α_norm[i] · η_bare^{-α_norm[i]}
    other_subs = Dict{Any,Any}()
    for (i, xi) in enumerate(indep_vars)
        i == pivot_idx && continue
        iszero(alpha_norm[i]) && continue
        other_subs[xi] = pivot^alpha_norm[i] * η_bare^(-alpha_norm[i])
    end
    expr = S(expr, other_subs)

    # 5. Substitute pivot = 1 (PDE is scaling-homogeneous so pivot drops out).
    expr2 = S(expr, Dict(pivot => Symbolics.Num(1)))

    # 6. Check no independent variables remain.
    remaining        = Num.(Symbolics.get_variables(expr2))
    still_has_indep  = any(r -> any(isequal(r), indep_vars), remaining)
    still_has_indep && return nothing

    # 7. Divide out the minimal power of η_bare if present.
    expr2 = _divide_min_eta_power(expr2, η_bare)

    return expr2
end


# ---------------------------------------------------------------------------
# Helper: simplify by dividing out the minimum η_bare power
# ---------------------------------------------------------------------------

"""
    _min_poly_degree(expr, var) -> Int

Compute the minimum power of `var` over all additive terms of `expr`.
For a monomial, this equals its degree. For a sum, it's the minimum over terms.
"""
function _min_poly_degree(expr, var)
    raw = Symbolics.unwrap(expr)
    if SymbolicUtils.iscall(raw) && SymbolicUtils.operation(raw) === (+)
        terms = SymbolicUtils.arguments(raw)
        degs = [try convert(Int, Symbolics.degree(Num(t), var)) catch; 0 end for t in terms]
        return minimum(degs)
    end
    return try convert(Int, Symbolics.degree(expr, var)) catch; 0 end
end

"""
    _termwise_divide(expr, η_bare, k) -> Num

Divide each additive term of `expr` by η_bare^k individually.
`simplify` fails to cancel η powers when the expression is a sum containing
Differential operators; working term-by-term avoids this.
"""
function _termwise_divide(expr, η_bare, k)
    k == 0 && return expr
    raw = Symbolics.unwrap(expr)
    if SymbolicUtils.iscall(raw) && SymbolicUtils.operation(raw) === (+)
        terms = SymbolicUtils.arguments(raw)
        new_terms = Num[simplify(Num(t) * η_bare^(-k); expand=true) for t in terms]
        return sum(new_terms)
    end
    return simplify(expr * η_bare^(-k); expand=true)
end

"""
    _divide_min_eta_power(expr, η_bare) -> Num

Divide out the minimum power of η_bare from the expression, simplifying
the ODE form.  Handles both plain sums and fraction (num/den) forms.
"""
function _divide_min_eta_power(expr, η_bare)
    raw = Symbolics.unwrap(expr)
    if SymbolicUtils.iscall(raw) && SymbolicUtils.operation(raw) === (/)
        # Fraction P(η)/η^d: factor η^k_min from numerator term-by-term
        num_raw, den_raw = SymbolicUtils.arguments(raw)
        k_min = _min_poly_degree(Num(num_raw), η_bare)
        k_min == 0 && return expr
        new_num = _termwise_divide(Num(num_raw), η_bare, k_min)
        new_den = simplify(Num(den_raw) * η_bare^(-k_min); expand=true)
        iszero(new_den - 1) && return new_num   # denominator cancelled completely
        return new_num / new_den
    end
    min_d = _min_poly_degree(expr, η_bare)
    min_d == 0 && return expr
    return _termwise_divide(expr, η_bare, min_d)
end


# ---------------------------------------------------------------------------
# Top-level: find_ode_dilation
# ---------------------------------------------------------------------------

"""
    find_ode_dilation(pde; indep_vars, dep_vars, verbose=false)

Dilation-method similarity solver.

Returns a Vector of result Dicts (same key structure as `find_ode`):
  "success", "PDE", "similarity_variable", "output_similarity", "PDE_similarity"
"""
function find_ode_dilation(pde; indep_vars::Vector, dep_vars::Vector, verbose::Bool=false)
    null_vecs, param_names = find_dilation_symmetry(pde, indep_vars, dep_vars)

    results = Dict{String,Any}[]

    for null_vec in null_vecs
        candidates = build_similarity_candidates(null_vec, param_names, indep_vars, dep_vars)
        for cand in candidates
            ode = reduce_to_ode(pde, indep_vars, dep_vars,
                                 cand.η_expr, cand.pivot, cand.pivot_idx,
                                 cand.gamma, cand.dep_idx, cand.alpha_norm)
            ode === nothing && continue

            dep = dep_vars[cand.dep_idx]
            ansatz_expr = iszero(cand.gamma) ? nothing : cand.pivot^cand.gamma

            result = Dict{String,Any}(
                "success"           => true,
                "PDE"               => ode,
                "similarity_variable" => cand.η_expr,
                "output_similarity"   => [dep] => ansatz_expr,
                "PDE_similarity"      => ode,
                "alpha_norm"          => cand.alpha_norm,
                "gamma"               => cand.gamma,
            )
            push!(results, result)
            verbose || break
        end
        verbose || (!isempty(results) && break)
    end

    return results
end


# ---------------------------------------------------------------------------
# String API: find_similarity_v2
# ---------------------------------------------------------------------------

"""
    find_similarity_v2(pde, bcs; parameters=[], verbose=false)

Like `find_similarity` but uses the dilation method internally.
Parses the same string syntax and returns the same result structure.
"""
function find_similarity_v2(
    pde::String,
    boundary_conditions::String;
    parameters::AbstractVector = [],
    verbose::Bool = false
)
    raw_bcs    = split(boundary_conditions, ';')
    bc_strings = [strip(s) for s in raw_bcs if !isempty(strip(s))]
    parsed_bcs = parse_boundary_condition.(bc_strings)

    input_vars = convert(Vector{Symbolics.Num},
        reduce(union, map(bc -> bc[2], parsed_bcs); init=Symbolics.Num[]))
    user_params = _normalize_parameters(parameters)
    bc_params   = convert(Vector{Symbolics.Num},
        reduce(union, map(bc -> bc[3], parsed_bcs); init=Symbolics.Num[]))
    all_params  = convert(Vector{Symbolics.Num}, union(user_params, bc_params))

    restrictions = map(bc -> bc[1], parsed_bcs)
    output_vars  = convert(Vector{Symbolics.Num},
        reduce(union, map(bc -> bc["function"], restrictions); init=Symbolics.Num[]))

    symbolic_pde = parse_pde(pde, input_vars, output_vars; parameters=all_params)

    results = find_ode_dilation(symbolic_pde;
                                 indep_vars=input_vars,
                                 dep_vars=output_vars,
                                 verbose=verbose)

    results = boundary_condition_similarity!(results, restrictions; input_vars=input_vars)

    if !verbose && !isempty(results)
        filtered = filter(p -> occursin("similarity", p[1]), results[1])
        return filtered
    end
    return results
end
