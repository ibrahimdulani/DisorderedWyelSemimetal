using Plots
using DataFrames
using CSV 
using LaTeXStrings

sigma=0.5
w_c = 0.24
nu=3
#phi=0


df = DataFrame(CSV.File("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/results.csv")) 

df = unique(df,[:L,:sig,:w])

ws = unique(df.w) 

plot(100^(1/nu).*(ws.-w_c) , filter([:L, :sig]=>(L,sig)-> L==100 && sig == sigma, df).alpha_q04,xlabel=L"L^{1/\nu}(w- w_c)",ylabel=L"\alpha_q", label="L=100",title="simga=$sigma, q=0.4, nu=$nu, Wc=$w_c", legend=:topleft, shape=:circle)

plot!(500^(1/nu).*(ws.-w_c) , filter([:L, :sig]=>(L,sig)-> L==500 && sig == sigma, df).alpha_q04, label="L=500", shape=:circle)

plot!(1000^(1/nu).*(ws.-w_c) , filter([:L, :sig]=>(L,sig)-> L==1000 && sig == sigma, df).alpha_q04, label="L=1000", shape=:circle)

plot!(5000^(1/nu).*(ws.-w_c) , filter([:L, :sig]=>(L,sig)-> L==5000 && sig == sigma, df).alpha_q04, label="L=5000", shape=:circle)

plot!(10000^(1/nu).*(ws.-w_c) , filter([:L, :sig]=>(L,sig)-> L==10000 && sig == sigma, df).alpha_q04, label="L=10000", shape=:circle)