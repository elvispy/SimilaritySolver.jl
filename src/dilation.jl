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
            # In newer Symbolics, op has an .order field (e.g. D^3)
            ord = try Int(op.order) catch; 1 end
            return _add_degrees(inner_deg, ScalingDegree(indep_map[dx_nm] => -Rational{Int}(ord)))
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
        unw_exp = Symbolics.unwrap(exp_raw)
        exp_val = if hasproperty(unw_exp, :val)
            unw_exp.val
        elseif unw_exp isa Real
            unw_exp
        else
            try Float64(unw_exp) catch; nothing end
        end
        
        (isnothing(exp_val) || !isfinite(Float64(exp_val))) && return _zero_degree()
        exp_rat = rationalize(Int, Float64(exp_val); tol=1e-9)
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
    pdes = (pde isa AbstractVector) ? pde : [pde]
    indep_map, dep_map = _make_maps(indep_vars, dep_vars)
    # Free symbols that are neither indep nor dep vars (physical parameters) get
    # zero scaling degree automatically in _scaling_degree_raw.
    # stable column ordering
    param_names = vcat([indep_map[k] for k in [Symbolics.unwrap(v).name for v in indep_vars]],
                       [dep_map[begin_dep_name(v)] for v in dep_vars])

    n_params = length(param_names)
    A_rows = Vector{Vector{Float64}}()

    for p in pdes
        raw   = Symbolics.unwrap(p)
        terms = _collect_add_terms(raw)
        degs  = [_scaling_degree_raw(t, indep_map, dep_map) for t in terms]
        n_terms = length(degs)
        n_terms <= 1 && continue

        ref = degs[1]
        for i in 2:n_terms
            diff = _add_degrees(degs[i], _neg_degree(ref))
            row = zeros(Float64, n_params)
            for (j, pname) in enumerate(param_names)
                row[j] = Float64(get(diff, pname, 0//1))
            end
            push!(A_rows, row)
        end
    end

    if isempty(A_rows)
        return zeros(Float64, 0, n_params), param_names
    end

    A = reduce(hcat, A_rows)' |> Matrix{Float64}
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
# Null-space combination search (Phase 2.5)
# ---------------------------------------------------------------------------

"""
    _expand_null_space_combinations(null_vecs; N=4) -> Vector{Vector{Rational{Int}}}

Given the basis vectors of a null space, return all unique integer linear combinations
`Σ cᵢ·vᵢ` with `|cᵢ| ≤ N`, normalized and de-duplicated.  The original basis vectors
are always included first.  For a 1-D null space returns the input unchanged.

This is needed because the physically relevant similarity variable often lies in the span
of the null-space basis but is not itself a basis vector (e.g. the Barenblatt-Pattle
solution for nonlinear diffusion, where the physical vector is `2v₁ - v₂`).
"""
function _expand_null_space_combinations(null_vecs::Vector{Vector{Rational{Int}}}; N::Int=4)
    k = length(null_vecs)
    k <= 1 && return null_vecs          # nothing to combine for 1-D null space

    seen   = Set{Vector{Rational{Int}}}()
    # Collect (L1_norm_of_coeffs, normalised_vector) for sorting
    pending = Vector{Tuple{Int, Vector{Rational{Int}}}}()

    # Helper: normalise a rational vector (first-nonzero positive, integer, primitive)
    function _normalise(v::Vector{Rational{Int}})
        all(iszero, v) && return nothing
        idx = findfirst(!iszero, v)
        v   = v[idx] < 0 ? -v : copy(v)
        den = lcm([x.den for x in v if !iszero(x)])
        vi  = [Rational{Int}(numerator(x * den), 1) for x in v]
        g   = gcd([abs(numerator(x)) for x in vi if !iszero(x)])
        return vi .// g
    end

    # Enumerate all integer-coefficient combinations (including basis vectors as unit coeffs)
    ranges = [(-N):N for _ in 1:k]
    for coeffs in Iterators.product(ranges...)
        all(iszero, coeffs) && continue
        l1 = sum(abs, coeffs)
        v_comb = sum(coeffs[i] * null_vecs[i] for i in 1:k)
        vn = _normalise(v_comb)
        vn === nothing && continue
        vn in seen && continue
        push!(seen, vn)
        push!(pending, (l1, vn))
    end

    # Sort by L1 norm of coefficients: simpler combinations first (parsimony order)
    sort!(pending; by = t -> t[1])
    return [t[2] for t in pending]
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
            exp_i = -(α_norm[i].den // α_norm[i].num)
            η_expr = η_expr * indep_vars[i]^exp_i
        end

        push!(results, (
            η_expr     = η_expr,
            pivot      = xp,
            pivot_idx  = pivot_idx,
            gamma_vals = γ_norm,
            alpha_norm = α_norm,
        ))
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
                           gamma::Rational{Int}, f_fn, η_expr,
                           orders::Vector{Int})
    pivot_var = indep_vars[pivot_idx]
    k_piv     = orders[pivot_idx]

    # Use a dummy function η_fn for the chain rule expansion
    η_fn_sym = eval(:(@variables η_dil(..)))[1]
    η_sym = η_fn_sym(indep_vars...)

    result = Symbolics.Num(0)
    for j in 0:k_piv
        ff = _ffact(gamma, j)
        iszero(ff) && continue
        u_j = ff * pivot_var^(gamma - j)
        
        # 1. Build derivative chain of f(η_sym)
        v = f_fn(η_sym)
        for i in eachindex(indep_vars)
            cnt = (i == pivot_idx) ? (k_piv - j) : orders[i]
            for _ in 1:cnt; v = Differential(indep_vars[i])(v); end
        end
        # 2. Expand: f(η_sym)_x = f'(η_sym) * η_sym_x
        v = expand_derivatives(v)
        # 3. Substitute η_sym -> η_expr and expand again: η_sym_x -> (x/t)_x = 1/t
        # Use diff2term to ensure derivatives are matched correctly
        v = Symbolics.substitute(Symbolics.diff2term(Symbolics.value(v)), 
                                 Dict(Symbolics.diff2term(Symbolics.value(η_sym)) => η_expr))
        v = expand_derivatives(v)
        # 4. Hide η_expr -> η_bare
        v = Symbolics.substitute(v, Dict(η_expr => η_bare))
        
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
    # Handle Vector of expressions
    if raw isa AbstractVector
        for r in raw
            append!(keys, _scan_dep_keys(r, dep_raw))
        end
        return keys
    end

    # Unwrap if Num
    raw_unw = Symbolics.unwrap(raw)
    dep_unw = Symbolics.unwrap(dep_raw)

    # Base case: is this node dep_raw itself?
    if isequal(raw_unw, dep_unw)
        push!(keys, raw_unw)
        return keys
    end
    
    SymbolicUtils.iscall(raw_unw) || return keys
    op   = SymbolicUtils.operation(raw_unw)
    args = SymbolicUtils.arguments(raw_unw)

    # Recursive check for all children (dep may appear nested inside products, etc.)
    for a in args
        append!(keys, _scan_dep_keys(a, dep_unw))
    end

    # If this is a differential chain, check if inner part is a dep-chain
    if op isa Differential
        inner_keys = _scan_dep_keys(args[1], dep_unw)
        if !isempty(inner_keys)
            push!(keys, raw_unw)   # this outer node is also a dep-chain key
        end
    end

    return unique(keys)
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
    _build_dep_subs(pdes, dep, indep_vars, pivot_idx, gamma, f_fn, η_sym)

Scan `pdes` for every derivative pattern applied to `dep`, then return a
substitution dictionary mapping each pattern to the correct Leibniz expansion
of  pivot^γ · f(η_sym).

This is order-agnostic: it works for PDEs of any derivative order without any
pre-specified maximum.
"""
function _build_dep_subs(pdes, dep, indep_vars::Vector, pivot_idx::Int,
                          gamma::Rational{Int}, f_fn, η_sym)
    dep_raw = Symbolics.unwrap(dep)

    # Collect dep-chain keys from all PDEs (scalar or vector input)
    all_raw_keys = if pdes isa AbstractVector
        mapreduce(p -> _scan_dep_keys(Symbolics.unwrap(p), dep_raw), vcat, pdes)
    else
        _scan_dep_keys(Symbolics.unwrap(pdes), dep_raw)
    end
    raw_keys = unique(all_raw_keys)

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
                        gamma_vals::Vector{Rational{Int}},
                        alpha_norm::Vector{Rational{Int}})
    # 1. Build substitutions for ALL dependent variables across ALL PDEs.
    # Expand derivatives first so that D(r*u*h) → u*h + r*Du*h + …
    pde_expanded = pde isa AbstractVector ? expand_derivatives.(pde) : expand_derivatives(pde)
    dep_subs = Dict{Any,Any}()
    f_fns = [eval(:(@variables $(Symbol("f_dil", j))(..)))[1] for j in eachindex(dep_vars)]
    for (j, dep) in enumerate(dep_vars)
        merge!(dep_subs, _build_dep_subs(pde_expanded, dep, indep_vars, pivot_idx, gamma_vals[j], f_fns[j], η_expr))
    end
    _sub_expand(e) = simplify(expand_derivatives(substitute(e, dep_subs)); expand=true)
    expr = pde_expanded isa AbstractVector ? _sub_expand.(pde_expanded) : _sub_expand(pde_expanded)

    # 2. Express every non-pivot variable through η_bare:
    #    xi = pivot^α_norm[i] · η_bare^{-α_norm[i]}
    other_subs = Dict{Any,Any}()
    for (i, xi) in enumerate(indep_vars)
        i == pivot_idx && continue
        iszero(alpha_norm[i]) && continue
        other_subs[xi] = pivot^alpha_norm[i] * η_bare^(-alpha_norm[i])
    end
    
    if expr isa AbstractVector
        expr = [simplify(substitute(e, other_subs); expand=true) for e in expr]
        expr = [simplify(substitute(e, Dict(pivot => Symbolics.Num(1))); expand=true) for e in expr]
        expr = [simplify(expand_derivatives(e); expand=true) for e in expr]
    else
        expr = simplify(substitute(expr, other_subs); expand=true)
        expr = simplify(substitute(expr, Dict(pivot => Symbolics.Num(1))); expand=true)
        expr = simplify(expand_derivatives(expr); expand=true)
    end

    # 4. Check no independent variables remain.
    expr2 = expr
    remaining = if expr2 isa AbstractVector
        Num.(mapreduce(e -> Symbolics.get_variables(e), union, expr2))
    else
        Num.(Symbolics.get_variables(expr2))
    end
    still_has_indep  = any(r -> any(isequal(r), indep_vars), remaining)
    still_has_indep && return nothing

    # 5. Divide out the minimal power of η_bare if present.
    if expr2 isa AbstractVector
        return [_divide_min_eta_power(e, η_bare) for e in expr2]
    else
        return _divide_min_eta_power(expr2, η_bare)
    end
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

Dilation-method similarity solver.  Accepts a single PDE (`Num`) or a coupled
system (`Vector{Num}`); in the latter case the shared dilation symmetry is found
by stacking constraint rows from all PDEs into one invariance matrix.

Returns a `Vector{Dict{String,Any}}`. Each Dict has keys:
- `"success"` — `true`
- `"PDE"` / `"PDE_similarity"` — reduced ODE (`Num`) for scalar input, or
  `Vector{Num}` (one ODE per input PDE) for system input
- `"similarity_variable"` — the similarity variable η
- `"output_similarity"` — `dep_vars => ansatz_exprs`
- `"alpha_vec"` — `Vector{Rational}` of normalized indep-variable scaling exponents
- `"gamma_vals"` — `Vector{Rational}` of scaling exponents, one per dep_var
- `"gamma"` — scalar shorthand for `gamma_vals[1]` (backward compat; single dep_var)
"""
function find_ode_dilation(pde; indep_vars::Vector, dep_vars::Vector, verbose::Bool=false)
    null_vecs, param_names = find_dilation_symmetry(pde, indep_vars, dep_vars)

    # Expand with integer linear combinations so physically relevant similarity
    # variables not aligned with a basis vector are also found (Phase 2.5).
    all_vecs = _expand_null_space_combinations(null_vecs; N=4)

    results = Dict{String,Any}[]

    for null_vec in all_vecs
        candidates = build_similarity_candidates(null_vec, param_names, indep_vars, dep_vars)
        for cand in candidates
            ode = reduce_to_ode(pde, indep_vars, dep_vars,
                                 cand.η_expr, cand.pivot, cand.pivot_idx,
                                 cand.gamma_vals, cand.alpha_norm)
            ode === nothing && continue

            ansatz_exprs = Num[iszero(g) ? Num(1) : cand.pivot^g for g in cand.gamma_vals]

            result = Dict{String,Any}(
                "success"             => true,
                "PDE"                 => ode,
                "similarity_variable" => cand.η_expr,
                "output_similarity"   => dep_vars => ansatz_exprs,
                "PDE_similarity"      => ode,
                "alpha_vec"           => cand.alpha_norm,
                "gamma_vals"          => cand.gamma_vals,
                # backward-compat scalar key (first dep_var); for systems use "gamma_vals"
                "gamma"               => cand.gamma_vals[1],
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
    find_similarity_v2(pde, boundary_conditions; parameters=[], verbose=false)

String-API wrapper around [`find_ode_dilation`](@ref) that mirrors the interface
of `find_similarity` (the original heuristic solver) but uses the exact
dilation/scaling-symmetry method instead.

# Arguments
- `pde::String`: PDE in string notation, e.g. `"du/dt + 6*u*du/dx + d3u/d3x = 0"`.
  Derivatives are written as `d<n><var>/d<n><wrt>` where `<n>` is the order
  (omit for first order). The right-hand side after `=` is subtracted automatically.
- `boundary_conditions::String`: Semicolon-separated boundary conditions, e.g.
  `"u(x=Inf, t) = 0; u(x, t=0) = 1"`. Each entry fixes one or more independent
  variables and assigns a value to the dependent variable there.
- `parameters`: Optional list of parameter names (String or `Symbolics.Num`) that
  appear in the PDE but are not independent or dependent variables (e.g. `["ν", "U∞"]`).
- `verbose::Bool`: If `false` (default), returns a filtered `Dict` of the first
  similarity result containing only keys with `"similarity"` in their name
  (`"similarity_variable"`, `"PDE_similarity"`, `"output_similarity"`).
  If `true`, returns the full `Vector{Dict}` from `find_ode_dilation`.

# Returns
- `verbose=false`: `Dict{String,Any}` with keys
  `"similarity_variable"`, `"PDE_similarity"`, `"output_similarity"`.
- `verbose=true`: `Vector{Dict{String,Any}}` — one entry per null-space candidate,
  each with all keys from `find_ode_dilation` (including `"PDE"`, `"gamma"`,
  `"substitutions"`, `"similarity_variable"`, etc.).
- Empty `Dict` / empty `Vector` if no similarity solution exists.

# Example
```julia
using SimilaritySolver
result = find_similarity_v2(
    "du/dt - d2u/d2x = 0",
    "u(x=Inf, t) = 0; u(x, t=0) = 1"
)
println(result["similarity_variable"])  # η = x * t^(-1//2)
println(result["PDE_similarity"])       # reduced heat ODE: f'' + (1/2)η f' = 0
```
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

    # Extract dep vars from the PDE string, NOT from BCs.
    # parse_boundary_condition misparses derivative BCs like "dψ/dy(x,y=0)=0"
    # as a new function "dy", which would poison dep_vars with a spurious variable,
    # widen the null space, and crash reduce_to_ode.
    dep_names = unique([m[1] for m in eachmatch(r"d\d?(\w+)(?:\(.*?\))?/", pde)])
    indep_strs = string.(input_vars)
    param_strs = string.(all_params)
    dep_names  = filter(n -> !(n in indep_strs) && !(n in param_strs), dep_names)
    output_vars = convert(Vector{Symbolics.Num}, [
        eval(:(@variables $(Symbol(n))(..)))[1](input_vars...)
        for n in dep_names
    ])

    symbolic_pde = parse_pde(pde, input_vars, output_vars; parameters=all_params)

    results = find_ode_dilation(symbolic_pde;
                                 indep_vars=input_vars,
                                 dep_vars=output_vars,
                                 verbose=verbose)

    # NOTE: boundary_condition_similarity! (from similaritySolve.jl) expects the old
    # find_ode result format (output_similarity as a Pair, substitutions as a Dict of
    # similarity var → expr). The new find_ode_dilation format is incompatible.
    # BC transformation is skipped here for now; tracked as a roadmap item.
    # results = boundary_condition_similarity!(results, restrictions; input_vars=input_vars)

    if !verbose && !isempty(results)
        filtered = filter(p -> occursin("similarity", p[1]), results[1])
        return filtered
    end
    return results
end
