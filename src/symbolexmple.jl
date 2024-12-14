using Symbolics
using SymbolicUtils: arguments, addterms
# @variables R a A Γ ρ U Q Q1

# Q = 0.5A*R*a*Γ*ρ*(U^2)
# Q_string::String = repr(Q)

# # how to convert Q_string back into a symbolic expression of type Num?
# # this is not doing what I want to
# Q1 = eval(Meta.parse(Q_string))

# tree = Symbolics.unwrap(Q)
# tree1 = Symbolics.unwrap(Q1)
# @assert istree(tree)

# # this assertion fails
# @assert istree(tree1)


# @variables x y u(x, y) η(x, y) v(η) n::Int
# ddx = Differential(x); ddy = Differential(y)
# # Define the PDE
# expr = ddx(u) #- ddy(ddy(u))
# # Substitute u(x, y) with v(η)
# f_eta = substitute(expr, u => v)
# # Substitute η(x, y) = y / x back into f_eta
# f_eta_substituted = substitute(f_eta, η => y / x)
# println(f_eta_substituted)

using Symbolics

@variables x y v(x, y) η(x, y) f(η) n::Int m::Int
ddx = Differential(x); ddy = Differential(y)

# Define η = y / x^m
η_expr = y * x^m

# Define the PDE
v = y^n * f
expr =  expand_derivatives(ddy(ddy(v))) # - ddx(v)
println(expand_derivatives(expr))
# Expand the chain rule manually
dη_dx = Differential(x)(η_expr)  # Compute ∂η/∂x
dη_dy = Differential(y)(η_expr)  # Compute ∂η/∂y

expr = expand_derivatives(substitute(expr, Differential(x)(η) =>  dη_dx))
expr = expand_derivatives(substitute(expr, Differential(y)(η) =>  dη_dy))
println("----------------")
println(expr)
println("---")
# Step 1: Substitute η(x, y) => η to make it independent
expr = expand_derivatives(substitute(expr, η => :η))
println(expr)
println("---")

# Step 2: Replace y in terms of η and x
expr = substitute(expr, y => η * x^(-m))

# Expand to clean up derivatives
expr = expand_derivatives(expr)
expr = simplify(expr);
expr = substitute(expr, Dict([n => 1, m => -1/2]));
println(simplify(expr==0; expand=true));

if istree(expr) && operation(expr) == +
    terms = arguments(expr)
else
    terms = [expr]
end

# Print the decomposed terms
println("Decomposed terms: ", terms)

any(isequal(x), Symbolics.get_variables(simplify(expr*x^(.5);expand=true)))