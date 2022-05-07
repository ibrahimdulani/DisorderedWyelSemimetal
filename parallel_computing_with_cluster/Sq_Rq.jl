using Random
using Distributions
using FFTW
using LinearAlgebra 
using Arpack 


L=parse(Int64, ARGS[1])
sigma=parse(Float64, ARGS[2]); 
W=parse(Float64, ARGS[3])

theta=0.5;  #the twist

Nbox=25;   #number of boxes (which should be selected such that L/Nbox must be integer)
q = [0.1, 0.4];        #exponent in generalized inverse participation ratio 
Rq=zeros(length(q))
Sq=zeros(length(q))
NW = 20;    # number of disorder strengths
ws=0.05; 

#Random.seed!(234) #setting the seed 


tk(n,sig,ll,th)=-sign(cos(2.0*pi*n/ll+th/ll))*abs(cos(2.0*pi*n/ll+th/ll))^sig; #k-space Hamiltonian 


Vx = W.*randn(L) 
Vk = fft(Vx)/L 

Ham = zeros(Complex{Float64},L,L) 

for i in 1:L
    for j in 1:L
        if i==j
            Ham[i,j]=tk(i,sigma,L,theta) #+ Vkk/L
        elseif i<j
            Ham[i,j]=Vk[j-i+1]
        else
            Ham[i,j]=conj(Vk[i-j+1])
         end
     end
 end 

 ev0, evec0 = eigs(Ham, nev=1, which=:SM) 

 Nr = sum(abs.(evec0).^2); 

 mu=zeros(Nbox) 

 for k in 1:Nbox
    mu[k]=sum(abs.(evec0[((k-1)*div(L,Nbox))+1:k*div(L,Nbox)]).^2)/Nr;
 end

 for i in 1:length(q)
    Rq[i]=sum(mu.^q[i]);
    Sq[i]=sum((mu.^q[i]).*(log.(mu)));
 end

println("#L: ", L, ", sigma: ", sigma, ", W: ", W)

print("#")
 for i in 1:length(q)
    print("S",q[i],",","R",q[i],",")
 end
 println()
 for i in 1:length(q)
    print(Sq[i],",",Rq[i],",")
 end
 println()
