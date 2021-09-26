using LinearAlgebra

println("k_local(Elastic Modulus,thickness,height,length,number of nodes,row one, row two, column one, column two)")

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