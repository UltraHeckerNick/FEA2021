using LinearAlgebra

# k_local(Elastic Modulus,thickness,height,length,number of nodes,row/column one, row/column two)

#creates local stiffness matrix and maps to global
function k_local(k,l,n,rc_1,rc_2)
    k_global = zeros(n,n)
    k_global[rc_1,rc_1] = k
    k_global[rc_1,rc_2] = -k
    k_global[rc_2,rc_1] = -k
    k_global[rc_2,rc_2] = k
    return k_global
end

# constants from problem statement
k1 = 40
k2 = 30
k3 = 60
k4 = 50
# combining all globally mapped matrices
kG = k_local(k1,l,5,1,2)+
     k_local(k2,l,5,2,4)+
     k_local(k3,l,5,3,4)+
     k_local(k4,l,5,4,5)

# Adjusts Global Matrix for fixed-end reactions
function adj_fix(kG,x,n)
    default = zeros(Float64,1,n)
    default[x] = 1
    kG[x,:] = default
    kG[:,x] = transpose(default)
end
adj_fix(kG,1,5)
adj_fix(kG,3,5)

# Force Vector
f = Float64[0 0 0 100 60]
F = transpose(f)

# Displacement vector
u = kG\F
println("Nodal displacement vector (mm):")
display(u)


