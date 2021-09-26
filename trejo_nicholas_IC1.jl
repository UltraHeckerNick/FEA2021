using LinearAlgebra
using Symbolics

# Problem One

A = Matrix{Int8}([0 5 0; 8 3 7; 9 -2 9])
B = Matrix{Int8}([4 6 -2; 7 2 3; 1 3 -4])
C = Vector{Int8}([-1, 2, 5])

println("Result of A*B: ")
display(A*B)
println("")

println("Result of A*C: ")
display(A*C)
println("")

println("Identity Matrix of Matrix A: ")
identityA = Matrix{Int8}(I,3,3)
display(identityA)
println("")

println("Result of I*A: ")
display(identityA*A)
println("")

println("Result of A*I: ")
display(A*identityA)
println("")

println("Matrix A")
display(A)
println("")

#Problem Two

A = Matrix{Int8}([5  7  2  0  3  5;
     3  8 -3 -5  0  8;
     1  4  0  7 15  9;
     0 10  5 12  3 -1;
     2 -5  9  2 18 -10])

B = Matrix{Int8}([ 2 10  0;
      8  7  5;
     -5  2 -4;
      4  8 13;
      3 12  0;
      1  5  7])

println("Result C = A*B: ")
display(A*B)
println("")

println("Transpose of Matrix A: ")
display(transpose(A))
println("")

# Problem Three

A = Matrix{Int8}([1 5 0; 8 3 7; 6 -2 9])

println("Determinant of Matrix A: ")
display(det(A))
println("")

println("Determinant of the Transpose of Matrix A: ")
display(det(transpose(A)))
println("")

println("Inverse of Matrix A: ")
display(inv(A))
println("")

# Problem Four
@variables k1
A = [0 0 0 0;
     0 k1 0 -k1;
     0 0 0 0;
     0 -k1 0 k1;]
println("Local Matrix A: ")
display(A)
println("")

G = [0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0;
     0 0 0 k1 0 0 0 -k1 0 0;
     0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0;
     0 0 0 -k1 0 0 0 k1 0 0;
     0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0;]
println("Mapped global matrix, G: ")
display(G)









