# Plan: JOSS Publication Roadmap for SimilaritySolver.jl

## Current State (as of 2026-03-10)

### What Works ✓
- `find_ode_dilation` — exact rational dilation symmetry + PDE→ODE reduction
- **System-PDE support**: `find_ode_dilation(pdes::Vector{Num}; ...)` stacks constraints from all PDEs
- **EDS 2025 validation**: IM-regime coupled system (mass/momentum/surfactant) — invariance matrix rank 2, EDS null vector verified, 3 coupled ODEs produced
- **Parsimony ordering**: results sorted by L1 norm of null-space combination coefficients (basis vectors first)
- **`gamma_vals`**: Vector{Rational} for multi-dep-var systems; backward-compat `"gamma"` key kept
- Heat equation, KdV, Burgers, 3-var heat, Barenblatt (nonlinear diffusion) all verified
- **Phase 4 Scaffold**: Documenter.jl set up with index, api, algorithm, comparison, and limitations pages.
- 284/284 tests pass (12 testsets)
- 3 exported functions: `find_similarity`, `find_similarity_v2`, `find_ode_dilation`

### Known Gaps
- Documentation lacks "Architectural Decisions" and "Community Guidelines" required by 2026 JOSS standards.
- `boundary_condition_similarity!` incompatible with new dict format (BC transform skipped in `find_similarity_v2`)

---

## Literature Survey Summary (arxiv 2022–2025)

**Field state**: Boundary layer, KdV, and nonlinear diffusion similarity solutions are still
being computed **by hand** in 2025. No open-source Python or Julia tool does end-to-end
automated similarity reduction. Maple (GeM/PDEtools) and Mathematica (MathLie) dominate,
but they require manual user input for symmetry identification.

**Closest competition**: Druzhkov & Cheviakov's "Invariant Reduction" series
(arXiv:2412.02965, 2501.09313, 2507.08213 — 3 papers in 2024–2025) implements
systematic PDE→ODE reduction in Maple. It does NOT auto-detect dilation symmetries;
users still provide the symmetry group manually. Julia/Python binding: none.

**The unique position of SimilaritySolver**:
*Only tool that takes a symbolic PDE, automatically detects dilation symmetry via
exact linear algebra, and returns the reduced ODE — in Julia, in one function call.*

---

## JOSS Roadmap (Revised Order)

**Phase order**: Phase 4 (Documentation & JOSS Compliance) → Phase 3 (JOSS Paper) → Phase 5 (Ecosystem)

---

### Phase 4 — Documentation & JOSS Compliance (est. 3–5 days)
*Documenter.jl setup and JOSS "Scope and Significance" requirements.*

1. **Architectural Decisions Doc** — Document the "why" behind the design: exact rational arithmetic vs floats, the Leibniz workaround for Symbolics.jl issues, and the decision to extract `dep_vars` from the PDE string.
   - [ ] 4.1 Create `docs/src/architecture.md`
   - [ ] 4.2 Document trade-offs (Performance vs. Mathematical Rigor)

2. **Community & Contribution Guidelines** — JOSS requires "pathways for community contribution."
   - [ ] 4.3 Create `docs/src/contributing.md`
   - [ ] 4.4 Define bug reporting process and feature request workflow.
   - [ ] 4.5 Add a `COMMUNITY.md` or update `index.md` with clear "How to help" section.

3. **Refine Statement of Need** — Explicitly articulate the "problem framing" and why SimilaritySolver is a "substantial contribution" to the research community.
   - [ ] 4.6 Update `docs/src/comparison.md` with citable research impact.

4. **Installation & Verification (The "Colleague Test")** — Ensure installation and test execution are foolproof.
   - [ ] 4.7 Review and update `docs/src/getting_started.md`.

---

### Phase 3 — JOSS Paper (est. 2–3 days)
*Write the paper once Phase 4 JOSS compliance is verified.*

5. **Write `paper.md`** — JOSS format. Sections:
    - Summary
    - Statement of Need
    - Algorithm (point to docs)
    - Functionality (heat, KdV, Blasius examples)
    - Comparison (table)
    - References

6. **Add `paper.bib`** — BibTeX entries for Maple GeM, SymPy, the Cheviakov series, EDS 2025, etc.

---

### Phase 5 — Ecosystem Integration (post-JOSS)
*High-impact integrations that attract users from adjacent Julia ecosystems.*

7. **ModelingToolkit.jl `PDESystem` input** — optional extension dispatch.
8. **NeuralPDE.jl tutorial notebook** — MTK PDE → `find_ode_dilation` → ODE solve → PINN compare.
9. **Julia Registry registration** — `]add SimilaritySolver`.

---

## Master Checklist

### Phase 4 — Documentation & JOSS Compliance

#### 4.1 Architectural Decisions & Trade-offs
- [ ] 4.1a Create `docs/src/architecture.md`
- [ ] 4.1b Document the rationale for exact rational arithmetic (avoiding numerical drift in scaling)
- [ ] 4.1c Document the Leibniz workaround as a principled design choice for high-order PDEs
- [ ] 4.1d Explain the `dep_vars` parsing strategy (PDE extraction vs. BC extraction)
- [ ] 4.1e Include "Human Contributions": how we frame the PDE reduction problem

#### 4.2 Community & Contribution Guidelines
- [ ] 4.2a Create `docs/src/contributing.md`
- [ ] 4.2b Add "Reporting Bugs" section (without referencing internal bug trackers)
- [ ] 4.2c Add "Development Setup" section for new contributors
- [ ] 4.2d Add "Code of Conduct" or link to standard Julia CoC

#### 4.3 Statement of Need & Impact
- [ ] 4.3a Update `docs/src/comparison.md` to emphasize research utility and "citability"
- [ ] 4.3b Contrast with Cheviakov's GeM (Maple) regarding automation level

---

### Phase 3 — JOSS Paper

- [ ] 3.1 Create `paper.md` with JOSS YAML front matter
- [ ] 3.2 Write Statement of Need (2026 JOSS style)
- [ ] 3.3 Create `paper.bib` with all relevant citations
- [ ] 3.4 Run JOSS paper validator

---

## Verification After Each Phase

- **Phase 4**: `julia --project=docs docs/make.jl` builds with all new sections; manual review for JOSS compliance.
- **Phase 3**: `joss paper` validator passes.
