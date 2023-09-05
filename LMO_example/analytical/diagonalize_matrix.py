from sympy import Matrix, I, simplify
from sympy.abc import A, B

# Define the matrix
M = Matrix([
    [A, 0, 0, -B],
    [0, A, -B, 0],
    [0, B, -A, 0],
    [B, 0, 0, -A]
])

# Compute the diagonal matrix and the diagonalizing matrix
P, D = M.diagonalize()

# Normalize the columns of P
for col in range(P.cols):
    P[:, col] = P[:, col].normalized()

# Ensure the matrices are unitary: P*P.H should be the identity, where P.H is the conjugate transpose
is_unitary = simplify(P * P.H) == Matrix.eye(4)
assert is_unitary, "Matrix P is not unitary!"

# Compute the inverse of P
P_inv = P.inv()

# Display results
print("Diagonalized matrix (D):")
print(D)
print("\nDiagonalizing matrix (P):")
print(P)
print("\nInverse of the diagonalizing matrix (P^-1):")
print(P_inv)
