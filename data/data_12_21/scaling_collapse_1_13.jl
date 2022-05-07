using DataFrames 
using CSV
using Plots
using LaTeXStrings
using Statistics 
using Random 

df = DataFrame(CSV.File("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/results.csv")) 

#df = df[completecases(df),:]

miss_or_nan_indices = [] 

# miss_or_nan_indices = append!(miss_or_nan_indices,findall(isnan,df.alpha_q0))
# miss_or_nan_indices = append!(miss_or_nan_indices,findall(isnan,df.alpha_q04))
# miss_or_nan_indices = append!(miss_or_nan_indices,findall(isnan,df.alpha_q05))
# miss_or_nan_indices = append!(miss_or_nan_indices,findall(isnan,df.alpha_q06))
# miss_or_nan_indices = append!(miss_or_nan_indices,findall(isnan,df.alpha_q08))
# miss_or_nan_indices = append!(miss_or_nan_indices,findall(isnan,df.alpha_q15))

# df = delete!(df,miss_or_nan_indices)

W(x,y,d,x1,y1,d1,x2,y2,d2)=((y-(((x2-x)*y1 - (x1-x)*y2)/(x2-x1)))/(sqrt(d^2 + (((x2-x)*d1)/(x2-x1))^2 +(((x1-x)*d2)/(x2-x1))^2)))^2 

#W(x,y,d,x1,y1,d1,x2,y2,d2)=(1/(sqrt(d^2 + (((x2-x)*d1)/(x2-x1))^2 +(((x1-x)*d2)/(x2-x1))^2)))^2 

sigs = sort(unique(df.sig)) 

Smin = zeros(Float64,length(sigs),2)

for s in 1:length(sigs)

    sigma=sigs[s]

    ww=filter([:L,:sig,:w]=>(L,sig,w)->L!=100 && sig==sigma && w!=0.01 && w!=0.1,df).w
    ll=filter([:L,:sig,:w]=>(L,sig,w)->L!=100 && sig==sigma && w!=0.01 && w!=0.1,df).L

    wcs=[minimum(ww)+0.1:0.005:maximum(ww)-0.1;]

    nus=[-10:0.1:10;]
    deleteat!(nus,findall(x->x==0.0,nus)) #delete 0.0 from the list

    S=zeros(Float64,length(wcs),length(nus))

    for i in 1:length(wcs)

        wc=wcs[i]

        for j in 1:length(nus)
            
            nu=nus[j]

            X=(ww.-wc).*(ll.^(1/nu))
            Y=(filter([:L,:sig,:w]=>(L,sig,w)->L!=100 && sig==sigma && w!=0.01 && w!=0.1,df).alpha_q06).*ll
            Dy=filter([:L,:sig,:w]=>(L,sig,w)->L!=100 && sig==sigma && w!=0.01 && w!=0.1,df).err06

            wsum = 0 
            for i in 2:length(X)-1
                wsum = wsum + W(X[i],Y[i],Dy[i],X[i-1],Y[i-1],Dy[i-1],X[i+1],Y[i+1],Dy[i+1])
            end

            S[i,j] = wsum/(length(X)-2)

        end

    end

    Smin[s,1]= wcs[argmin(S)[1]]
    Smin[s,2]= nus[argmin(S)[2]]

end 

plot(sigs,Smin[:,1], xlabel=L"\sigma", ylabel=L"W_c", title="q=0.6", label="", shape=:circle)