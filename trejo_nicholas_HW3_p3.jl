using LinearAlgebra

# k_local(stiffness coefficient,number of nodes,row/column one, row/column two, ==1 adjusts for fixed node)

#creates local stiffness matrix and maps to global
function k_local(k,i,n1,n2,fixed=0)
      k_global = zeros(i,i)
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
      display(k_global)
      return k_global
  end

# constants from problem statement
k1 = 5
k2 = 8
k3 = 5
k4 = 20
k5 = 10
k6 = 20
# combining all globally mapped matrices
kG = k_local(k1,5,1,2,1)+
     k_local(k2,5,2,3)+
     k_local(k3,5,2,3)+
     k_local(k4,5,2,4)+
     k_local(k5,5,3,4)+
     k_local(k6,5,5,4,1)

# Force Vector
f = Float64[0 10 0 10 0]
F = transpose(f)

# Displacement vector
u = kG\F
println("Nodal displacement vector (inches):")
display(u)

# Strain Function
ϵ(x2,x1,l) = (x2 - x1)/l

#Reaction Forces vector
Fapp(k,x2,x) = k*(x2-x)
FappM = [Fapp(k1,u[2],u[1])
         Fapp(k2,u[3],u[2])
         Fapp(k3,u[3],u[2])
         Fapp(k4,u[4],u[2])
         Fapp(k5,u[4],u[3])
         Fapp(k6,u[5],u[4])]
println("Reaction Forces vector (lbs):")
display(FappM)

