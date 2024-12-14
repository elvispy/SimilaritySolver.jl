using Symbolics
using SymbolicUtils

const SymbolicType = Union{Num, SymbolicUtils.BasicSymbolic}

"""
    is_output(var::Num, vars::Vector{Num})::Bool

Determines if the variable `var` is dependent on any other variable in the `vars` list.

# Arguments:
- `var::Num`: A symbolic variable to check.
- `vars::Vector{Num}`: A vector of symbolic variables.

# Returns:
- `true` if `var` depends on any other variable in `vars`.
- `false` otherwise.

The dependency is determined by checking if the derivative of `var` with respect to another variable is non-zero.
"""
function is_output(var::T, vars::Vector{T})::Bool where {T<:SymbolicType}
    for other_var in vars
        if other_var !== var  # Avoid self-comparison
            # If the derivative is non-zero, `var` depends on `other_var`.
            if (expand_derivatives(Differential(other_var)(var)) != 0) !== false
                return true
            end
        end
    end
    return false
end

"""
    is_input(var::Num, vars::Vector{Num})::Bool

Determines if the variable `var` acts as an input to any other variable in `vars`.

# Arguments:
- `var::Num`: A symbolic variable to check.
- `vars::Vector{Num}`: A vector of symbolic variables.

# Returns:
- `true` if `var` acts as an input to any variable in `vars`.
- `false` otherwise.

The input role is identified if the derivative of any variable in `vars` with respect to `var` is non-zero.
"""
function is_input(var::T, vars::Vector{T})::Bool where {T<:SymbolicType}
    for other_var in vars
        if other_var !== var  # Avoid self-comparison
            # If the derivative is non-zero, `var` is an input to `other_var`.
            if (expand_derivatives(Differential(var)(other_var)) != 0) !== false
                return true
            end
        end
    end
    return false
end

"""
    decomposeVars(vars::Vector{Num})

Splits a list of symbolic variables into inputs and outputs.

# Arguments:
- `vars::Vector{Num}`: A vector of symbolic variables.

# Returns:
- `(inputs, outputs)`: A tuple where:
  - `inputs`: Variables identified as independent inputs.
  - `outputs`: Variables identified as dependent outputs.

Variables are classified based on their dependencies:
- Outputs depend on other variables.
- Inputs influence the outputs.
"""
function decomposeVars(vars::Vector{T}) where {T<:SymbolicType}
    outputs = filter(var -> is_output(var, vars), vars)  # Identify dependent variables
    inputs  = filter(var -> is_input(var, outputs), vars)  # Identify independent variables
    return inputs, outputs
end

"""
    guess_powers(expression::Num; variables::Vector{Num}, values::Vector{AbstractVector})

(Placeholder function) Attempts to guess powers of variables that simplify the given symbolic expression.

# Arguments:
- `expression::Num`: The symbolic expression to analyze.
- `variables::Vector{Num}`: A list of variables whose powers are to be guessed.
- `values::Vector{AbstractVector}`: A list of possible values for the powers.

# Notes:
- The implementation for this function is currently not provided.
"""
function guess_powers(expression::T; variables::Vector{T}, values::Vector{AbstractVector}) where {T<:SymbolicType}
    # Implementation is not provided in the original script.
    powers_guess = collect(Iterators.product(values...));
    result = Dict{String, Any}()
    result["success"] = false;
    for guess = powers_guess
        dict_substitutions = Dict(zip(variables, guess));
        expression_guess = simplify(substitute(expression, dict_substitutions); expand=true);
        d = Symbolics.degree(expression_guess, sym=x);
        PDE = Symbolics.get_variables(simplify(expression*x^(-d);expand=true));
        if !any(isequal(x), PDE)
            result["success"] = true;
            result["PDE"] = PDE;
            result["substitutions"] = dict_substitutions;
            return results
        end
    end
    return results
end

"""
    findODE(symbolicPDE::Num)::Num

Finds a similarity solution to a given symbolic partial differential equation (PDE).

# Arguments:
- `symbolicPDE::Num`: A symbolic representation of the PDE.

# Returns:
- A simplified ordinary differential equation (ODE) if a similarity solution is found.
- `nothing` if no similarity solution is identified.

# Steps:
1. Decomposes the variables in the PDE into inputs and outputs.
2. Defines similarity variables (e.g., η = y * x^m).
3. Substitutes the similarity variables into the PDE.
4. Attempts to guess powers `n` and `m` that simplify the resulting equation.
"""
function findODE(symbolicPDE::T)::T where {T<:SymbolicType}
    vars = Num.(Symbolics.get_variables(symbolicPDE))  # Extract variables from the PDE
    inputs, outputs = decomposeVars(vars)        # Split variables into inputs and outputs
    parameters = setdiff(vars, inputs, outputs)  # Parameters are those not classified as inputs or outputs

    @assert length(inputs) == 2 "Only two independent variables are supported"
    @assert length(outputs) == 1 "Only one dependent variable is supported"

    # Define auxiliary similarity variables
    x, y = inputs
    @variables n::Int m::Int η(x, y) f(η)
    v = expand_derivatives(outputs[1])
    η_expr = y * x^m  # Define similarity variable η

    # Compute derivatives of η with respect to x and y
    dη_dx = Differential(x)(η_expr)
    dη_dy = Differential(y)(η_expr)

    # Substitute and simplify the PDE
    expr = expand_derivatives(substitute(symbolicPDE, v => y^n * f))
    expr = expand_derivatives(substitute(expr, Differential(x)(η) => dη_dx))
    expr = expand_derivatives(substitute(expr, Differential(y)(η) => dη_dy))
    expr = expand_derivatives(substitute(expr, η => :η))
    expr = expand_derivatives(substitute(expr, y => η * x^(-m)))

    # Attempt to guess powers for the solution
    results = guess_powers(expr; variables=[n, m], values=[[-1, 0, 1], -.5:.5:.5])
    
    if results["success"] == true
        return results.PDE  # Return the simplified ODE
    end

    println("Tests inconclusive. No similarity solution was found with guess η=$η_expr")
    return nothing
end