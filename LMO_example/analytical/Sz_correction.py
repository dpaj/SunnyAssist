from sympy import Matrix, sqrt, symbols

# Define the symbols
A, B = symbols('A B')

# Define the matrix X_q
X_q = Matrix([
    [(A - sqrt(A**2 - B**2)) / B, 0, (A + sqrt(A**2 - B**2)) / B, 0],
    [0, (A - sqrt(A**2 - B**2)) / B, 0, (A + sqrt(A**2 - B**2)) / B],
    [0, 1, 0, 1],
    [1, 0, 1, 0]
])

# Define the matrix X_q^†
X_q_dagger = Matrix([
    [-B / (2 * sqrt(A**2 - B**2)), 0, 0, (A * sqrt(A**2 - B**2) + A**2 - B**2) / (2 * (A**2 - B**2))],
    [0, -B / (2 * sqrt(A**2 - B**2)), (A * sqrt(A**2 - B**2) + A**2 - B**2) / (2 * (A**2 - B**2)), 0],
    [B / (2 * sqrt(A**2 - B**2)), 0, 0, (-A * sqrt(A**2 - B**2) + A**2 - B**2) / (2 * (A**2 - B**2))],
    [0, B / (2 * sqrt(A**2 - B**2)), (-A * sqrt(A**2 - B**2) + A**2 - B**2) / (2 * (A**2 - B**2)), 0]
])

print("X_q:")
print(X_q)
print("\nX_q^†:")
print(X_q_dagger)

A = Matrix([
    [1, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 0]
])

alpha = simplify(X_q_dagger * A * X_q)

print('transformed')
print(alpha)