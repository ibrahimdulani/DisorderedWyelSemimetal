using DataFrames 
using CSV
using Plots
using LaTeXStrings
using Statistics 
using Random 

df = DataFrame(CSV.File("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/results.csv")) 

W(x,y,d,x1,y1,d1,x2,y2,d2)=(y-(((x2-x)*y1 - (x1-x)*y2)/(x2-x1)))^2/(d^2 + (((x2-x)*d1)/(x2-x1))^2 +(((x1-x)*d2)/(x2-x1))^2) 

#sigs = sort(unique(df.sig)) 

sigma=0.1

ww=filter([:L,:sig]=>(L,sig)->L!=100 && sig==sigma,df).w 
ll=filter([:L,:sig]=>(L,sig)->L!=100 && sig==sigma,df).L

Y=(filter([:L,:sig]=>(L,sig)->L!=100 && sig==sigma,df).alpha_q04) 

# ll = 1000

# ww = filter([:L,:sig]=>(L,sig)->L==ll && sig==sigma,df).w
# Y = filter([:L,:sig]=>(L,sig)->L==ll && sig==sigma,df).alpha_q04

Dy = Y.-mean(Y)

wcs=[0.01:0.01:0.5;] 

nus=[0.05:0.5:10;]

S = zeros(length(wcs),length(nus))

for n in 1:length(wcs)

    for m in 1:length(nus)

        X = (ww.-wcs[n]).*(ll.^nus[m])

        ws = zeros(length(X)-2)

        for i in 2:length(X)-1
            ws[i-1] = W(X[i],Y[i],Dy[i],X[i-1],Y[i-1],Dy[i-1],X[i+1],Y[i+1],Dy[i+1])
        end 

        S[n,m] = mean(ws)

    end

end


# Smin = zeros(Float64,length(sigs),2)

# for s in 1:length(sigs)

#     sigma=sigs[s]

#     ww=filter([:L,:sig]=>(L,sig)->L!=100 && sig==sigma,df).w 
#     ll=filter([:L,:sig]=>(L,sig)->L!=100 && sig==sigma,df).L

#     wcs=[minimum(ww):0.05:maximum(ww);]

#     nus=[0:0.1:5;]
#     deleteat!(nus,findall(x->x==0.0,nus)) #delete 0.0 from the listS

#     S=zeros(Float64,length(wcs),length(nus))

#     for k in 1:length(wcs)

#         wc=wcs[k]

#         for j in 1:length(nus)
            
#             nu=nus[j]

#             X=(ww.-wc).*(ll.^(1/nu))
#             Y=(filter([:L,:sig]=>(L,sig)->L!=100 && sig==sigma,df).alpha_q04)#.*ll
#             Dy=Y.-mean(filter(!isnan,Y)) 

#             wsum = 0 
#             for i in 2:length(X)-1
#                 wsum = wsum + W(X[i],Y[i],Dy[i],X[i-1],Y[i-1],Dy[i-1],X[i+1],Y[i+1],Dy[i+1])
#             end

#             S[k,j] = wsum/(length(X)-2)

#         end

#     end

#     Smin[s,1]= wcs[argmin(S)[1]]
#     Smin[s,2]= nus[argmin(S)[2]]

# end 

# plot(sigs,Smin[:,1], xlabel=L"\sigma", ylabel=L"W_c", title="q=0.4", label="", shape=:circle)