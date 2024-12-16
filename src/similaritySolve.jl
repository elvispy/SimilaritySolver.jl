using Symbolics
using SymbolicUtils
using Logging

const SymbolicType = Union{Num, SymbolicUtils.BasicSymbolic, Nothing}

iszerosym(sym) = (simplify(sym; expand=true) != 0) === false

"""
    is_output(var::T, vars::Vector{T})::Bool where {T<:SymbolicType}

Check if the symbolic variable `var` depends on any variable in the list `vars`.

# Arguments
- `var::T`: A symbolic variable of type `SymbolicType` to check.
- `vars::Vector{T}`: A vector of symbolic variables to compare dependencies against.

# Returns
- `true` if `var` depends on any variable in `vars` (excluding itself).
- `false` otherwise.

# Details
The function computes the partial derivative of `var` with respect to each variable in `vars`
(using `expand_derivatives` and `Differential`). If any derivative is **non-zero** and 
the variable is not `var` itself, the function returns `true`.

# Example
```julia
using Symbolics

@variables x y z
vars = [x, y]

f = x^2 + y + z
is_output(f, vars)  # true

g = x^2 + z
is_output(g, [y])   # false
is_output(g, [x])   # true
"""
function is_output(var::T, vars::Vector{T})::Bool where {T<:SymbolicType}
    for other_var in vars
        if other_var !== var  # Avoid self-comparison
            # If the derivative is non-zero, `var` depends on `other_var`.
            if !iszerosym(expand_derivatives(Differential(other_var)(var)))
                return true
            end
        end
    end
    return false
end


"""
    is_input(var::T, vars::Vector{T})::Bool where {T<:SymbolicType}

Check if the symbolic variable `var` acts as an input to any variable in the list `vars`.

# Arguments
- `var::T`: A symbolic variable of type `SymbolicType` to check.
- `vars::Vector{T}`: A vector of symbolic variables to check against.

# Returns
- `true` if `var` acts as an input to any variable in `vars` (excluding itself).
- `false` otherwise.

# Details
The function computes the partial derivative of each variable in `vars` with respect to `var`
(using `expand_derivatives` and `Differential`). If any derivative is **non-zero**, the function 
returns `true`.

# Example
```julia
using Symbolics

@variables x y z
vars = [x, y]

f = x^2 + y
is_input(x, [f])  # true
is_input(z, [f])  # false
"""
function is_input(var::T, vars::Vector{T})::Bool where {T<:SymbolicType}
    for other_var in vars
        if other_var !== var  # Avoid self-comparison
            # If the derivative is non-zero, `var` is an input to `other_var`.
            if !iszero(expand_derivatives(Differential(var)(other_var))) != false
                return true
            end
        end
    end
    return false
end

"""
    decomposeVars(vars::Vector{T})::Tuple{Vector{T}, Vector{T}} where {T<:SymbolicType}

Split a list of symbolic variables into inputs and outputs.

# Arguments
- `vars::Vector{T}`: A vector of symbolic variables of type `SymbolicType`.

# Returns
- `(inputs, outputs)`: A tuple containing:
- `inputs`: Variables that act as independent inputs (influencing other variables).
- `outputs`: Variables that depend on other variables.

# Details
The function classifies variables as:
- **Outputs**: Variables that depend on other variables in `vars` (using `is_output`).
- **Inputs**: Variables that influence the outputs (using `is_input`).

# Example
```julia
using Symbolics

@variables x y z
vars = [x, y, z]

# Decompose variables into inputs and outputs
inputs, outputs = decomposeVars(vars)

println("Inputs: ", inputs)   # Variables influencing others
println("Outputs: ", outputs) # Variables depending on others
"""
function decomposeVars(vars::Vector{T}) where {T<:SymbolicType}
    outputs = filter(var -> is_output(var, vars), vars)  # Identify dependent variables
    inputs  = filter(var -> is_input(var, outputs), vars)  # Identify independent variables
    return inputs, outputs
end

"""
    guess_powers(expression::Num; variables::Vector{Num}, values::Vector{AbstractVector})::Dict{String, Any}

Attempts to guess powers of variables that simplify the given symbolic expression into a solvable PDE.

# Arguments
- `expression::Num`: The symbolic expression to analyze.
- `variables::Vector{Num}`: A list of independent variables whose powers are to be guessed.
- `values::Vector{AbstractVector}`: A list of possible values for the powers.

# Returns
A `Dict{String, Any}` with the following keys:
- `"success"::Bool`: Indicates whether the guess was successful.
- `"PDE"::Num`: The simplified Partial Differential Equation (PDE) if successful.
- `"substitutions"::Dict{Any, Any}`: A dictionary containing the guessed power substitutions if successful.
# Notes
- The function iterates through all possible combinations of the provided values for the powers of the variables.
- It substitutes these guesses into the expression and checks if the resulting PDE has degree zero in the primary variable.
- If a valid substitution is found, it simplifies the PDE and returns the result.
"""
function guess_powers(expression::T; indep_vars::Vector{T}, powers::Vector{T}, values) where {T<:SymbolicType}
    powers_guess = collect(Iterators.product(values...));
    result = Dict{String, Any}()
    result["success"] = false;
    x, y, η = indep_vars;
    for guess = powers_guess
        dict_substitutions = Dict(zip(powers, guess));
        @debug "Trying $dict_substitutions"
        expression_guess = simplify(substitute(expression, dict_substitutions); expand=true);
        d = convert(Float64, Symbolics.degree(expression_guess, x));
        PDE = simplify(expression_guess*x^(-d);expand=true);
        @debug "Result: $PDE"
        if maximum(abs.(Symbolics.degree.([PDE, substitute(PDE, x => 1/x)], x))) == 0.0
            if any(var -> isequal(x)(var) || isequal(y)(var), Symbolics.get_variables(PDE))
                try
                    @variables aux
                    PDE2 = substitute(PDE, η => aux)
                    PDE2 = substitute(PDE2, Dict(x => 1.0, y=> 1.0))
                    PDE = simplify(substitute(PDE2, aux => η); expand=true);
                    PDE = simplify(substitute(PDE, Dict(i => convert(Float64, i) for i in -5:5)); expand=true);
                catch
                    @warn "Degree is zero but couldn't remove independent variable"
                end
            end
            result["success"] = true;
            result["PDE"] = PDE;
            result["substitutions"] = dict_substitutions;
            @info "Got an ODE with $dict_substitutions"
            return result
        end
    end
    return result
end

function trySimilarity(output_expr, η_expr, symbolicPDE, similarity_var, inputs, outputs, powers);
    x, y = inputs;
    η = similarity_var;
    v = expand_derivatives(outputs[1]);

    #η_expr = y * x^m  # Define similarity variable η
    @info "Trying general similarity guess $(outputs[1]) = $(substitute(output_expr, η=> η_expr))..."
    # Compute derivatives of η with respect to x and y
    dη_dx = Differential(x)(η_expr)
    dη_dy = Differential(y)(η_expr)

    # Substitute and simplify the PDE
    SES(input, dict) = simplify(expand_derivatives(substitute(input, dict)); expand=true);
    expr = SES(symbolicPDE, v => output_expr);
    expr = SES(expr, Dict(Differential(x)(η) => dη_dx, Differential(y)(η) => dη_dy));
    expr = SES(expr, η => :η);
    xd = Symbolics.degree(η_expr, x); yd = Symbolics.degree(η_expr, y);
    expr = SES(expr, y => η^(1/yd) / x^(xd/yd));
    expr = SES(expr, :η => η);
    
    #@info "Done."
    # Attempt to guess powers for the solution
    @info "Choosing the right values for the exponents..."
    unique_fractions = n -> sort(collect(-numerator / denominator for denominator in 1:n for numerator in 1:abs(denominator)-1 if gcd(numerator, abs(denominator)) == 1))
    results = guess_powers(expr; indep_vars=[x, y, η], powers=powers, values=[[0.0, 1.0, -2/3], unique_fractions(6)])
    return results
end


"""
    find_ode(symbolicPDE::Num)::Num

Finds a similarity solution to a given symbolic partial differential equation (PDE).

# Arguments
- `symbolicPDE::Num`: A symbolic representation of the partial differential equation (PDE).

# Returns
- A simplified ordinary differential equation (ODE) if a similarity solution is found.
- `nothing` if no similarity solution is identified.

# Steps
1. Decomposes the variables in the PDE into inputs (independent variables) and outputs (dependent variables).
2. Defines similarity variables (e.g., η = y * x^m) based on the inputs.
3. Substitutes the similarity variables into the PDE.
4. Attempts to guess powers `n` and `m` for the similarity transformation that simplify the resulting equation.
5. If a valid transformation is found, returns the simplified ODE.
6. If no solution is found, returns `nothing`.

# Notes
- The function assumes that the PDE involves two independent variables and one dependent variable.
- The function uses a trial-and-error approach to guess the appropriate powers for the similarity variable transformation, with a range of values considered for the exponents.

# Example
```julia
using Symbolics

@variables x y u(x, y)
symbolicPDE = Differential(x)(u) - Differential(y)(Differential(y)(u))
vars = [x, y, u];
result = find_ode(symbolicPDE; vars=vars)
if result !== nothing
    println("Similarity solution found: ", result)
else
    println("No similarity solution found.")
end

"""
function find_ode(symbolicPDE::T; vars::Vector{T}=nothing, log::Bool=true) where {T<:SymbolicType}
    if vars === nothing
        vars = Num.(Symbolics.get_variables(symbolicPDE))  # Extract variables from the PDE
    end
    log ? nothing : Logging.disable_logging(Logging.Info)
    inputs, outputs = decomposeVars(vars)        # Split variables into inputs and outputs
    parameters = setdiff(vars, Num.(inputs), outputs)  # Parameters are those not classified as inputs or outputs

    @assert length(inputs)  == 2 "Only two independent variables are supported"
    @assert length(outputs) <= 2 "At most two dependent variables are supported"

    # Define auxiliary similarity variables
    x, y = inputs; 
    @variables n::Float64 m::Float64 η($x, $y) f(η)
    η_exprs = [y * x^m, x * y^m];
    output_expr = [y^n * f, x^n * f];
    result = nothing; results = [];
    for (η_exp, out_exp) = Iterators.product(η_exprs, output_expr)
        result = trySimilarity(out_exp, η_exp, symbolicPDE, η, inputs, outputs, [n, m]);
        if result["success"] == true
            result["out_exp"] = outputs => out_exp
            result["similarity_variable"] = η => substitute(η_exp, result["substitutions"]);
            result["value"] = Symbolics.diff2term(Symbolics.value(Symbolics.substitute(result["PDE"]/η^Symbolics.degree(result["PDE"], η), η => Symbolics.variable("η"))));

            @info "Got similarity with f=$out_exp, η=$η_exp"
            push!(results, result)
            #break
        else
            @info "Similarity unsuccessful."
        end
        @info "-------------"
    end
    @info "Done."

    if length(results) == 0
        @info "Tests inconclusive. No similarity solution was found with guesses $(String(Symbol(outputs[1])))=$output_expr"
    end

    Logging.disable_logging(Logging.Debug)

    return results
end


"""
    parse_boundary_condition(condition::String)::Vector{Any}

Parses a boundary condition string into its function, restrictions, and value.

# Arguments
- `condition::String`: A boundary condition string of the format `"f(var1=val1, var2) = value"`. 
  The function name, variable restrictions, and value must follow this structure.

# Returns
A `Vector` containing:
1. A `Dict` with the following keys:
   - `"function"`: A symbolic function with its variables.
   - `"restriction"`: A `Dict{Symbol, Union{Float64, Num, Nothing}}` containing the variable restrictions (e.g., `x=0`).
   - `"value"`: The numeric or symbolic value associated with the boundary condition.
2. A list of symbolic variables (`Vector{Num}`) corresponding to the input variables.

# Description
- Uses a regular expression to parse the boundary condition into its components:
  - Extracts the function name, variable restrictions, and the condition's value.
- Parses the restriction string into a dictionary of variable-value pairs.
- Creates symbolic variables for the function and its inputs using `Symbolics`.
- Handles numeric or symbolic values for the condition.

# Notes
- Restrictions may omit values, such as `"t"` in `f(x=0, t) = 1.0`. These are parsed as `nothing`.
- If the value cannot be parsed as a `Float64`, it is interpreted as a symbolic variable.

# Examples
```julia
using Symbolics

# Example 1: Parsing a boundary condition with numeric values
bc = "f(x=0, t) = 1.0"
result = parse_boundary_condition(bc)

println(result[1]["function"])      # f(x, t)
println(result[1]["restriction"])   # Dict(:x => 0.0, :t => nothing)
println(result[1]["value"])         # 1.0

# Example 2: Parsing a boundary condition with symbolic value
bc2 = "g(y=1) = alpha"
result2 = parse_boundary_condition(bc2)

println(result2[1]["function"])      # g(y)
println(result2[1]["restriction"])   # Dict(:y => 1.0)
println(result2[1]["value"])         # alpha
"""
function parse_boundary_condition(condition)
    # Use a regular expression to extract the function, restriction, and value
    regex = r"(?<function>\w+)\((?<restriction>[^\)]+)\)\s*=\s*(?<value>[\s\-\w\.]+)"
    
    # Extract the match for the boundary condition
    match_result = match(regex, condition)
    value = nothing
    # If a match is found
    if match_result !== nothing
        # Extract components using named capturing groups
        restriction_str = match_result["restriction"]
        try
            value = parse(Float64, match_result["value"])
        catch
            value = Symbolics.variable(strip(match_result["value"]))
        end
        
        # Convert the restriction string into pairs (e.g., "x=0, t" => [(:x, 0), (:t, 0)])
        restriction = Dict{Union{Symbol, Num}, Union{Float64, Num, Nothing}}()
        input_vars = [];
        for r in split(restriction_str, ",")
            var, val = strip.(split(r, "=")) ∪ [nothing]
            varNum = eval(Meta.parse("@variables $var"))[1]
            restriction[varNum] = isnothing(val) ? val : parse(Float64, val)
            push!(input_vars, Symbolics.variable(var));
        end
        eval(Meta.parse("@variables "*join(input_vars, " ")))
        function_value = eval(Meta.parse("@variables $(strip(match_result["function"]))($(join(input_vars, ", ")))"))[1];
        #function_value = Symbolics.variable(, T=Symbolics.FnType)(input_vars...);
        
        return [Dict("function" => function_value, "restriction" => restriction, 
                       "value" => value),  input_vars]
    else
        throw(ArgumentError("Invalid boundary condition format: $condition"))
    end
end

"""
    nested_apply(input_string, f, start)

Applies a function `f` iteratively to variables extracted from the input string `input_string`. The function parses variable-multiplier pairs in the format `d<multiplier><variable>`, where `multiplier` is an integer (optional, defaults to 1) and `variable` is a string. The function then applies `f` to the variables in reverse order, starting with the value `start`.

### Arguments:
- `input_string::String`: A string containing variable-multiplier pairs in the format `d<multiplier><variable>`, where the multiplier is optional (defaults to 1) and the variable is a string.
- `f::Function`: A function that takes two arguments, the current result and the next variable, and returns the updated result.
- `start`: The initial value to start the iteration with, which will be updated by successive applications of `f`.

### Returns:
- The result of applying `f` to the variables in reverse order, starting from `start`.
"""
function nested_apply(input_string, f, start)
    # Step 1: Extract all variable-multiplier pairs
    pattern = r"d(\d*)([^d\s]+)"  # Matches d<multiplier><variable>
    matches = eachmatch(pattern, input_string)
    
    if isempty(matches)
        error("Invalid input format: $input_string")
    end
    
    # Step 2: Parse variables and their multipliers
    vars = []
    for m in matches
        multiplier = isempty(m[1]) ? 1 : parse(Int, m[1])  # Default multiplier to 1 if empty
        variable = m[2]
        for _ in 1:multiplier
            push!(vars, variable)  # Repeat variable 'multiplier' times
        end
    end
    
    # Step 3: Build the nested function call
    result = start  # Start with 0 as the innermost second argument
    for var in reverse(vars)  # Process variables in reverse order
        result = f(result, var)
    end
    
    return result
end

"""
    parse_pde(expr::String, input_vars, output_vars)::Num

Parse a partial differential equation (PDE) expression and convert differential terms into symbolic representations.

# Arguments
- `expr::String`: The PDE expression in string format, potentially containing differential terms (e.g., `d<multiplier><variable>`).
- `input_vars::Vector{Symbol}`: A vector of input variables for the PDE, used to define symbolic representations.
- `output_vars::Vector{Symbol}`: A vector of output variables for the PDE, used to define symbolic representations.

# Returns
- A numerical or symbolic result after parsing and evaluating the PDE expression.

# Details
This function replaces differential terms in `expr` with their symbolic derivatives using the `nested_apply` function. The input and output variables are declared symbolically before evaluating the modified PDE expression.

# Examples
```julia-repl
expr = "d2x/dy = x * y"
input_vars = [:x, :y]
output_vars = [:u]
parse_pde(expr, input_vars, output_vars)
"""
function parse_pde(expr::String, input_vars, output_vars)::Num
    # Replace '=' with negative form to prepare the expression
    expr = replace(expr, r"=(.+)" => s"- ( \1 )")
    
    # Apply differential transformation to captured terms
    apply_der(s, s_apply) = "(Differential($s_apply))($s)"
    for rgm = eachmatch(r"d\d?(\w+)(\(.*?\))?\/(\w+)", expr)
        equivalent_expr = nested_apply(rgm[3], apply_der, rgm[1])
        expr = replace(expr, rgm.match => equivalent_expr)
    end
    
    # Declare input and output variables
    eval(Meta.parse("@variables " * join(input_vars, " ")))
    eval(Meta.parse("@variables " * join(output_vars, " ")))

    # Evaluate the final expression
    return eval(Meta.parse(expr))
end

function boundary_condition_similarity!(results, restrictions; input_vars)
    if length(results) > 0
        for ii = eachindex(results)
            var = results[ii]
            x, y = input_vars; 
            @variables n::Float64 m::Float64 η($x, $y) f(η)
            output_function = Symbolics.substitute(var["out_exp"][2], var["similarity_variable"]);
            output_function = simplify(Symbolics.substitute(output_function, var["substitutions"]));
            BC_similarity = [];
            for restriction = restrictions
                push!(BC_similarity, Symbolics.substitute(output_function, filter(var -> !isnothing(var[2]), restriction["restriction"])) == restriction["value"]);
            end
            results[ii]["BC_similarity"] = BC_similarity;
        end
    else
        @warn "No similarity was found"
    end
    return results
end

function find_similarity(pde::String, boundary_conditions::String)

    # Split the boundary conditions by semicolons
    boundary_list = String.(split(boundary_conditions, ";"))

    boundary_conditions = parse_boundary_condition.(boundary_list)
    input_vars = convert(Vector{Num}, union(map(var -> var[2], boundary_conditions)...));
    restrictions = map(var -> var[1], boundary_conditions);
    output_vars = convert(Vector{Num}, union(map(var -> var["function"], restrictions)...));
    #println(iszero(Differential(input_vars[1])(output_vars[1])))
    #output_vars_dict = Dict(map(var -> chop("$(var["function"])", tail=1) => var["function"], restrictions));
    
    symbolicPDE = parse_pde(pde, input_vars, output_vars);
    results = find_ode(symbolicPDE; vars=input_vars ∪ output_vars);
    results = boundary_condition_similarity!(results, restrictions; input_vars=input_vars);
    # info = Dict("input_vars" => input_vars, 
    #            "output_vars" => output_vars, 
    #           "restrictions" => restrictions,
    #           "analysis_result" => results)
    return results
end

# Test the function with an example
#parsed_conditions = find_similarity("dw/dt = d2w/d2r", "w(r=Inf, t) = 0; w(r, t=0) = 1.0")
#println(parsed_conditions)

# @variables t r u(r, t)# η(x, y) f(η)
# symbolicPDE = Differential(t)(u) - Differential(r)(Differential(r)(u))
# vars = [t, r, u];
# result = find_ode(symbolicPDE; vars=vars)
# @variables x t u(t, x)# η(x, y) f(η)
# Dt = Differential(t); Dx = Differential(x)
# symbolicPDE = Dt(u) + 6*u*Dx(u) + Dx(Dx(Dx(u)))
# vars = [t, x, u];
# result = find_ode(symbolicPDE; vars=vars)