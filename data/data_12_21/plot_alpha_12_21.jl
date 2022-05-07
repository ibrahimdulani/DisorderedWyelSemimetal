using Plots
using DataFrames
using CSV 
using LaTeXStrings

df = DataFrame(CSV.File("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/results.csv"))

s=0.5
#st=string(s)

#q1=15
q=0.6

df = unique(df,[:L,:sig,:w])

plot(filter([:L, :sig]=>(L,sig)-> L==100 && sig == s, df).w[2:3],
    filter([:L, :sig]=>(L,sig)-> L==100 && sig == s, df).alpha_q06[2:3],label="L=100",
    xlabel="W",ylabel=L"\alpha_q",shape=:circle, title=L"\sigma=%$s, q=%$q",legend=:left)

plot!(filter([:L, :sig]=>(L,sig)-> L==500 && sig == s, df).w[2:3],
    filter([:L, :sig]=>(L,sig)-> L==500 && sig == s , df).alpha_q06[2:3],label="L=500",shape=:circle)

plot!(filter([:L, :sig]=>(L,sig)-> L==1000 && sig == s, df).w[2:3],
    filter([:L, :sig]=>(L,sig)-> L==1000 && sig == s , df).alpha_q06[2:3],label="L=1000",shape=:circle)

plot!(filter([:L, :sig]=>(L,sig)-> L==5000 && sig == s, df).w[2:3],
    filter([:L, :sig]=>(L,sig)-> L==5000 && sig == s, df).alpha_q06[2:3],label="L=5000",shape=:circle)

plot!(filter([:L, :sig]=>(L,sig)-> L==10000 && sig == s, df).w[2:3],
    filter([:L, :sig]=>(L,sig)-> L==10000 && sig == s, df).alpha_q06[2:3],label="L=10000",shape=:circle)
