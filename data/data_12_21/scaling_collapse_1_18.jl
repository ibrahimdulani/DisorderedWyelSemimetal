using DataFrames 
using CSV
using Plots
using LaTeXStrings
using Statistics 
using Random 
using NaNMath        #package for ignoring nan values

df = DataFrame(CSV.File("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/results.csv")) 

W(x,y,d,x1,y1,d1,x2,y2,d2)=((y-(((x2-x)*y1 - (x1-x)*y2)/(x2-x1)))/(sqrt(d^2 + (((x2-x)*d1)/(x2-x1))^2 +(((x1-x)*d2)/(x2-x1))^2)))^2 

sigs = sort(unique(df.sig)) 

Smin = zeros(Float64,length(sigs),2)

#sigma = 0.1 

for s in 1:length(sigs)

    sigma = sigs[s]
    
    ww=filter([:L,:sig]=>(L,sig)->L!=100 && L!=500 && sig==sigma,df).w 
    ll=filter([:L,:sig]=>(L,sig)->L!=100 && L!=500 && sig==sigma,df).L

    wcs = [minimum(ww)+0.05:0.005:maximum(ww)-0.05;] 
    nus = [-10:0.5:10;]
    deleteat!(nus,findall(x->x==0.0,nus)) #delete 0.0 from the list

    S=zeros(Float64,length(wcs),length(nus))

    
    for i in 1:length(wcs)
        
        wc = wcs[i] 

        for j in 1:length(nus)

            nu = nus[j]

            X=(ww.-wc).*(ll.^(1/nu))
            Y=filter([:L,:sig]=>(L,sig)->L!=100 && L!=500 && sig==sigma,df).alpha_q04
            Dy=filter([:L,:sig]=>(L,sig)->L!=100 && L!=500 && sig==sigma,df).err04

            # wc_inds = findall(x->x==wc,ww)

            # deleteat!(X,wc_inds) 
            # deleteat!(Y,wc_inds)
            # deleteat!(Dy,wc_inds)

            inds=sortperm(X)
            X=X[inds]
            Y=Y[inds]
            Dy=Dy[inds]

            wsum=zeros(length(X)-2)

            for k in 2:length(X)-1
                wsum[k-1] = W(X[k],Y[k],Dy[k],X[k-1],Y[k-1],Dy[k-1],X[k+1],Y[k+1],Dy[k+1])
            end 

            S[i,j]=NaNMath.mean(wsum)

        end

    end

    Smin[s,1]= wcs[argmin(S)[1]]
    Smin[s,2]= nus[argmin(S)[2]]

end

plot(sigs,Smin[:,1], xlabel=L"\sigma", ylabel=L"W_c", title="q=0.4", label="", shape=:circle)