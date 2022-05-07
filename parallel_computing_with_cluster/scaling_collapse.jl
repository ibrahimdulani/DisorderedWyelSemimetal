using Plots
using DataFrames
using CSV 
using LaTeXStrings

sigma=0.3
w_c = 0.34
nu=13 


df = DataFrame(CSV.File("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data/results.csv")) 

df = unique(df,[:L,:sig,:w])

ws = unique(df.w) 

plot(1000^(1/nu).*(ws.-w_c) ,filter([:L, :sig]=>(L,sig)-> L==1000 && sig == sigma, df).alpha_q1,xlabel=L"L^{1/\nu}(w- w_c)",ylabel=L"\alpha_q", label="L=1000",title=L"\sigma=0.3, q=0.2, \nu=13", legend=:topright, shape=:circle)

plot!(5000^(1/nu).*(ws.-w_c) ,filter([:L, :sig]=>(L,sig)-> L==5000 && sig == sigma, df).alpha_q1, label="L=5000", shape=:circle)

plot!(10000^(1/nu).*(ws.-w_c) ,filter([:L, :sig]=>(L,sig)-> L==10000 && sig == sigma, df).alpha_q1, label="L=10000", shape=:circle)