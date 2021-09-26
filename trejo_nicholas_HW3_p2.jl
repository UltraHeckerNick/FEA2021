using LinearAlgebra

# k_local(Elastic Modulus,thickness,height,length,number of nodes,row/column one,row/column two, ==1 adjusts for fixed node)

#creates local stiffness matrix and maps to global
function k_local(E,t,h,l,i,n1,n2,fixed=0)
    k_global = zeros(i,i)
    A = t*h
    k = ((E*A)/l)
    k_global[n1,n1] = k
    k_global[n1,n2] = -k
    k_global[n2,n1] = -k
    k_global[n2,n2] = k
    if fixed == 1
      adj = zeros(1,i)
      adj[n1] = 1
      k_global[n1,:] = adj 
      k_global[:,n1] = adj 
    end
    return k_global
end

# constants from problem statement
E = 28*10^3
t = 1/8
h16 = 4
h2345 = 2.5/4
l = 2
l2345 = 8

# combining all globally mapped matrices
kG = k_local(E,t,h16,l,4,1,2,1)+
     k_local(E,t,h2345,l2345,4,2,3)+
     k_local(E,t,h2345,l2345,4,2,3)+
     k_local(E,t,h2345,l2345,4,2,3)+
     k_local(E,t,h2345,l2345,4,2,3)+
     k_local(E,t,h16,l,4,3,4)

# Force Vector
f = Float64[0 0 0 500]
F = transpose(f)

# Displacement vector
u = kG\F
println("Nodal Displacement Vector (inches)")
display(u)

# Strain Function
ϵ(x2,x1,l) = (x2 - x1)/l

#Strain Vector
ϵM = [ϵ(u[2],u[1],l)
      ϵ(u[3],u[2],l2345) 
      ϵ(u[3],u[2],l2345) 
      ϵ(u[3],u[2],l2345) 
      ϵ(u[3],u[2],l2345)
      ϵ(u[4],u[3],l)]

# Stress Function
σ(E,ϵ) = E*ϵ

# Stress Vector
σM = [σ(E,ϵM[1])
      σ(E,ϵM[2])
      σ(E,ϵM[3])
      σ(E,ϵM[4])
      σ(E,ϵM[5])
      σ(E,ϵM[6])]
println("Element Stress Vector (psi)")
display(σM)

