using LinearAlgebra

#Givens from Problem Statement
#Elastic Modulus(psi),Area(in^2),length(in),Force(lb)
E = 1.9*10^6
A = 8
l = 3*12
l2 = 3*sqrt(2)*12
F = 500

# Stiffness Function
k(A,E,l) = (A*E)/l

# k_local(Elastic Modulus,Area,length,theta,number of nodes,
# row/column one,row/column two, ==1 left fixed/ ==2 right fixed)
function k_local(E,A,l,theta,i,n1,n2)
  p1 = (n1*2)-1
  p2 = n1*2
  p3 = (n2*2)-1
  p4 = n2*2
 
  k_global = zeros(2*i,2*i)
  kin = UniformScaling(k(A,E,l))
  s = (sind(theta))^2
	c = (cosd(theta))^2
	sc = sind(theta)*cosd(theta)
	k_quad = [c sc;
            sc s]
  k_global[p1:p2,p1:p2] = kin*k_quad
  k_global[p1:p2,p3:p4] = -kin*k_quad
  k_global[p3:p4,p1:p2] = -kin*k_quad
  k_global[p3:p4,p3:p4] = kin*k_quad
  
  return k_global
end

# Calculates element matrices
k_1 = k_local(E,A,l,0,5,1,2)
k_2 = k_local(E,A,l2,135,5,2,3)
k_3 = k_local(E,A,l,0,5,3,4)
k_4 = k_local(E,A,l,90,5,2,4)
k_5 = k_local(E,A,l2,45,5,2,5)
k_6 = k_local(E,A,l,0,5,4,5)

# Unadjusted Global
kG = k_1+k_2+k_3+k_4+k_5+k_6

# Constructs adjusted matrix; glb_adj(stiffness matrix, node number)
function glb_adj(k,n)
  nf = n*2
  n0 = nf-1
  k_global_adj = k[setdiff(1:end, (n0,nf)), setdiff(1:end, (n0,nf))] 
end
k_1_adj = glb_adj(k_1,1)
k_1_adj = glb_adj(k_1_adj,2)
k_2_adj = glb_adj(k_2,1)
k_2_adj = glb_adj(k_2_adj,2)
k_3_adj = glb_adj(k_3,1)
k_3_adj = glb_adj(k_3_adj,2)
k_4_adj = glb_adj(k_4,1)
k_4_adj = glb_adj(k_4_adj,2)
k_5_adj = glb_adj(k_5,1)
k_5_adj = glb_adj(k_5_adj,2)
k_6_adj = glb_adj(k_6,1)
k_6_adj = glb_adj(k_6_adj,2)


# Adjusted Global 
kG_adj = k_1_adj+k_2_adj+k_3_adj+k_4_adj+k_5_adj+k_6_adj

# Load Vector
p = Float64[0 0 0 -F 0 -F]
P = transpose(p)

# Displacement Vector
u = kG_adj\P
println("Displacement Vector (in)")
display(u)

# Reaction Vector
P_abs = kG*u_adj
println("Load Reaction Vector (lb):")
display(P_abs)

