using Symbolics
using LinearAlgebra
using Polynomials

@variables A B λ

# Define the matrix
M = [A 0 0 -B;
     0 A -B 0;
     0 B -A 0;
     B 0 0 -A]

# Calculate the characteristic polynomial
char_poly_expr = det(M - λ*I)

function extract_coefficients(expr, var)
    # If the polynomial is just a monomial or constant
    if typeof(expr) !== Symbolics.Add
        return [expr]  # Just return it as the coefficient
    end
    
    # We can obtain a degree by expanding the polynomial and counting the terms
    max_degree = maximum([term.exponent for term in Symbolics.arguments(expr) if hasproperty(term, :exponent)])
    coeffs = [0 for _ in 0:max_degree]

    for term in Symbolics.arguments(expr)
        if hasproperty(term, :exponent)
            coeff = term.coeff
            exp = term.exponent[1][2] # This assumes λ is the only variable with an exponent
            coeffs[exp+1] = coeff     # +1 because Julia is 1-based index
        else
            if term == var
                coeffs[2] = 1  # For terms like λ without an explicit coefficient
            else
                coeffs[1] = term
            end
        end
    end
    return coeffs
end

coeff_vector = extract_coefficients(char_poly_expr, λ)

# Display the coefficients
println("Characteristic Polynomial Coefficients:")
println(coeff_vector)

# Compute and display the eigenvalues
char_poly = Polynomial(reverse(coeff_vector))  # Construct polynomial
eigenvalues = Polynomials.roots(char_poly)  # Find its roots
println("Eigenvalues of the matrix M:")
println(eigenvalues)

