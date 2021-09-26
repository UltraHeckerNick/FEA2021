using LinearAlgebra

# k_local(Elastic Modulus,thickness,height,length,number of nodes,row one, row two, column one, column two)

#creates local stiffness matrix and maps to global
function k_local(E,t,h,l,n,r_1,r_2,c_1,c_2)
    k_global = zeros(n,n)
    A = t*h
    k = ((E*A)/l)
    k_global[r_1,c_1] = k
    k_global[r_1,c_2] = -k
    k_global[r_2,c_1] = -k
    k_global[r_2,c_2] = k
    return k_global
end

# constants from problem statement
E = 29*10^6
t = 1/16
h14 = 5
h23 = 2
l1 = 1
l23 = 4
l4 = 2

# combining all globally mapped matrices
kG = k_local(E,t,h14,l1,4,1,2,1,2)+
     k_local(E,t,h23,l23,4,2,3,2,3)+
     k_local(E,t,h23,l23,4,2,3,2,3)+
     k_local(E,t,h14,l4,4,3,4,3,4)
kG[1,:] = [1 0 0 0] #adjustment for Reaction from fixed end

# Force Vector
f = Float64[0 0 0 800]
F = transpose(f)

# Displacement vector
u = kG\F
println("Nodal displacement vector:")
display(u)

# Strain Function
ϵ(x2,x1,l) = (x2 - x1)/l

#Strain Vector
ϵM = [ϵ(u[2],u[1],l1)
      ϵ(u[3],u[2],l23) 
      ϵ(u[3],u[2],l23) 
      ϵ(u[4],u[3],l4) ]
println("Element strain vector")
display(ϵM)

# Stress Function
σ(E,ϵ) = E*ϵ

# Stress Vector
σM = [σ(E,ϵM[1])
      σ(E,ϵM[2])
      σ(E,ϵM[3])
      σ(E,ϵM[4])]
println("Element Stress Vector")
display(σM)

