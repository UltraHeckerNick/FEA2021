using LinearAlgebra

# k_local(Elastic Modulus,thickness,height,length,number of nodes,row one, row two, column one, column two)

#creates local stiffness matrix and maps to global
function k_local(k,l,n,r_1,r_2,c_1,c_2)
    k_global = zeros(n,n)
    k_global[r_1,c_1] = k
    k_global[r_1,c_2] = -k
    k_global[r_2,c_1] = -k
    k_global[r_2,c_2] = k
    return k_global
end

# constants from problem statement
E = 29*10^6
l = 2
k1 = 50
k2 = 60
k3 = 55
# combining all globally mapped matrices
kG = k_local(k1,l,4,1,2,1,2)+
     k_local(k2,l,4,2,3,2,3)+
     k_local(k3,l,4,3,4,3,4)
kG[4,:] = [0 0 0 1] #adjustment for Reaction from fixed end

# Force Vector
f = Float64[60 0 0 0]
F = transpose(f)

# Displacement vector
u = kG\F
println("Nodal displacement vector (inches):")
display(u)

# Strain Function
ϵ(x2,x1,l) = (x2 - x1)/l

#Strain Vector
ϵM = [ϵ(u[2],u[1],l)
      ϵ(u[3],u[2],l) 
      ϵ(u[4],u[3],l) ]
println("Element strain vector (inch/inch):")
display(ϵM)

# Stress Function
σ(E,ϵ) = E*ϵ

# Stress Vector
σM = [σ(E,ϵM[1])
      σ(E,ϵM[2])
      σ(E,ϵM[3])]
println("Element Stress vector (psi):")
display(σM)

