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
E = 29*10^6
l = 2
k1 = 25
k2 = 35
k3 = 45
k4 = 30
# combining all globally mapped matrices
kG = k_local(k1,l,5,1,2)+
     k_local(k2,l,5,2,3)+
     k_local(k3,l,5,3,4)+
     k_local(k4,l,5,3,5)
kG[4,:] = Float64[0 0 0 1 0] #adjustment for Reaction from fixed end
kG[5,:] = Float64[0 0 0 0 1] #adjustment for Reaction from fixed end
kG[:,4] = transpose(Float64[0 0 0 1 0]) #adjustment for Reaction from fixed end
kG[:,5] = transpose(Float64[0 0 0 0 1]) #adjustment for Reaction from fixed end

# Force Vector
f = Float64[75 0 0 0 0]
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
      ϵ(u[4],u[3],l) 
      ϵ(u[5],u[3],l)]
println("Element strain vector (inch/inch):")
display(ϵM)

# Stress Function
σ(E,ϵ) = E*ϵ

# Stress Vector
σM = [σ(E,ϵM[1])
      σ(E,ϵM[2])
      σ(E,ϵM[3])
      σ(E,ϵM[4])]
println("Element Stress vector (psi):")
display(σM)

#Reaction Forces vector
Fapp(k,x2,x) = k*(x2-x)
FappM = [Fapp(k12,u[2],u[1])
         Fapp(k12,u[3],u[2])
         Fapp(k3,u[4],u[3])
         Fapp(k4,u[5],u[3])]
println("Reaction Forces vector (lbs):")
display(FappM)

