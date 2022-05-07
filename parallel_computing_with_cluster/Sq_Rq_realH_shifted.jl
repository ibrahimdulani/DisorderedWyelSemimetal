using Random
using Distributions
using FFTW
using LinearAlgebra
using Arpack


sigma=parse(Float64,ARGS[1])
L=parse(Int64, ARGS[2])
W=parse(Float64, ARGS[3])

theta=0.5;  #the twist

Nbox=25;   #number of boxes (which should be selected such that L/Nbox must be integer)
q = [0.0, 0.4, 0.5, 0.6, 0.8, 1.5];        #exponent in generalized inverse participation ratio
Rq=zeros(length(q))
Sq=zeros(length(q))
#NW = 20;    # number of disorder strengths

#ws=0.05;

#Random.seed!(234) #setting the seed


tk(n,sig,ll,th)=-sign(cos(2.0*pi*n/ll+th/ll))*abs(cos(2.0*pi*n/ll+th/ll))^sig; #k-space Hamiltonian

Vx = W.*randn(L) 

#Vx = zeros(L)

Vxx = sum(Vx)/L

Hk = zeros(L)   #kinetic part of k-space Hamiltonian

for l in 1:L
    Hk[l] = tk(l,sigma,L,theta)
end 

Hf = ifft(Hk)

#Hf = ifftshift(Hf)

H = zeros(Complex{Float64},L,L); 

for i in 1:L
    for j in 1:L
        if i==j
            H[i,j] = Vx[i] - Vxx
        elseif i>j
            H[i,j] = Hf[i-j+1]
        else
            H[i,j] = conj(Hf[j-i+1]) 
        end
    end
end 

en0, evec0=eigs(H,nev=1,which=:SM);

Nrm = sqrt(sum(abs.(evec0).^2))

evec=evec0/Nrm

mu=zeros(Nbox)

for k in 1:Nbox
	mu[k]=sum(abs.(evec[((k-1)*div(L,Nbox))+1:k*div(L,Nbox)]).^2);
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
