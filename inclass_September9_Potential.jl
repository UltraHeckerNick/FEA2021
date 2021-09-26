using LinearAlgebra

# coefficients from PD of π functions
π = [70 -30 0;
     -30 140 -50;
     0 -50 50]
# Force Vector
f = Float64[0 100 60]
F = transpose(f)

# Displacement vector
u = π\F
println("Nodal displacement vector (mm):")
display(u)