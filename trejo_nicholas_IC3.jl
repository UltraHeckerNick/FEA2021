using LinearAlgebra

#Givens from Problem Statement
#Shear Modulus(psi),diameter(in),length(in),torque(lb*in)
G = 3.9*10^6
G2 = 4.0*10^6
d = 1.5
d2 = 1
l = 2.5*12
l2 = 1*12
l3 =2*12 
τ = -200*12 

# Polar Second Area Moment of a Cylinder (in^4)
J_cyl(d) = 0.5*pi*(d/2)^4
J = J_cyl(d)
J2 = J_cyl(d2)

# Torsional Stiffness Function
k(G,J,l) = (G*J)/l

# k_local(Elastic Modulus,thickness,height,length,number of nodes,
# row/column one,row/column two, ==1 left fixed/ ==2 right fixed)
function k_local(G,J,l,i,n1,n2,fixed=0)
    k_global = zeros(i,i)
    kin = k(G,J,l)
    k_global[n1,n1] = kin
    k_global[n1,n2] = -kin
    k_global[n2,n1] = -kin
    k_global[n2,n2] = kin
    if fixed == 1
      adj = zeros(1,i)
      adj[n1] = 1
      k_global[n1,:] = adj 
      k_global[:,n1] = adj 
    elseif fixed == 2
      adj = zeros(1,i)
      adj[n2] = 1
      k_global[n2,:] = adj 
      k_global[:,n2] = adj 
    end
    return k_global
end

# Global for angular displacement calculation
k_1 = k_local(G,J,l,4,1,2,1)
k_2 = k_local(G,J,l2,4,2,3)
k_3 = k_local(G2,J2,l3,4,3,4,2)
kG = k_1+k_2+k_3

#Global for Reaction calculation
k_11 = k_local(G,J,l,4,1,2)
k_22 = k_local(G,J,l2,4,2,3)
k_33 = k_local(G2,J2,l3,4,3,4)
kGG = k_11+k_22+k_33

# Torque Vector
t = Float64[0 τ 0 0]
T = transpose(t)

# Angular Displacement Vector
θ = kG\T
println("Nodal Angular Displacement Vector (radians)")
display(θ)

# Moment Reaction Vector
τ_app = kGG*θ
println("Moment Reaction Vector (lb*in):")
display(τ_app)
