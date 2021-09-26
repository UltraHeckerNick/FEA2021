using LinearAlgebra

#Givens from Problem Statement
#Elastic Modulus(psi),Area(in^2),length(in),Force(lb)
E = 1.9*10^6
A = 8
l = 3*12
l2 = 3*sqrt(2)*12
F = 500

# Torsional Stiffness Function
k(A,E,l) = (A*E)/l

# k_local(Elastic Modulus,Area,length,theta,number of nodes,
# row/column one,row/column two, ==1 left fixed/ ==2 right fixed)
function k_local(E,A,l,theta,i,n1,n2,fixed=0)
  #k_local = zeros(4,4)
  p1 = (n1*2)-1
  p2 = n1*2
  p3 = (n2*2)-1
  p4 = n2*2
 
  k_global = zeros(2*i,2*i)
  kin = k(A,E,l)
  s = (sind(theta))^2
	c = (cosd(theta))^2
	sc = sind(theta)*cosd(theta)
	k_quad = [c sc;
            sc s]
  k_global[p1:p2,p1:p2] = k*k_quad
  k_global[p1:p2,p3:p4] = -k*k_quad
  k_global[p3:p4,p1:p2] = -k*k_quad
  k_global[p3:p4,p3:p4] = k*k_quad
  # Need to remove rows and columns; check: push!, pop! 
  if fixed == 1
      adj = zeros(1,i)
      #adj[n1] = 1
      k_global[p1,:] = adj 
      k_global[:,p1] = adj 
      k_global[p2,:] = adj 
      k_global[:,p2] = adj 
    elseif fixed == 2
      adj = zeros(1,i)
      #adj[n2] = 1 
      k_global[p3,:] = adj 
      k_global[:,p3] = adj 
      k_global[p4,:] = adj 
      k_global[:,p4] = adj 
    end
    return k_global
end

# Unadjusted Global 
k_1 = k_local(E,A,l,0,5,1,2,1)
k_2 = k_local(E,A,l,135,5,2,3)
k_3 = k_local(E,A,l,0,5,3,4,1)
k_4 = k_local(E,A,l,90,5,2,4)
k_5 = k_local(E,A,l,45,5,2,5)
k_6 = k_local(E,A,l,0,5,4,5)
kG = k_1+k_2+k_3+k_4+k_5+k_6

#Adjusted Global 
kG_adj = kG

# Load Vector
p = Float64[0 0 0 -F 0 -F]
P = transpose(p)

# Angular Displacement Vector
u = kG_adj\P
println("Displacement Vector (in)")
display(u)

# Moment Reaction Vector
P_abs = kG*u_adj
println("Load Reaction Vector (lb):")
display(P_abs)
