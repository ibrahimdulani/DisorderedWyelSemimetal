using Printf
using Statistics
using CSV
using DataFrames 
using NaNMath        #package for ignoring nan values

L=5000
sigma=10
W=35
q = [0.0, 0.2, 0.5, 0.8, 1.2]
lambda=1/25


folder_name = @sprintf "/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_10_21/L%d/sigma%03d/W%03d/" L sigma W 

file_list = readdir(folder_name, join=true)

Sq = Array{Array{Float64}, 1}(undef, length(q))
Rq = Array{Array{Float64}, 1}(undef, length(q))

for j in 1:length(q)
	Sq[j] = Array{Float64, 1}(undef, length(file_list))
	Rq[j] = Array{Float64, 1}(undef, length(file_list))
end

for i in 1:length(file_list)
	str = readlines(file_list[i])
	if str==[] || length(str)<3
		skip
	else
		vals = split(str[3], ",")

		for j in 1:length(q)
			Sq[j][i] = parse(Float64, vals[2*j - 1])
			Rq[j][i] = parse(Float64, vals[2*j])
		end
	end
end

#print(size(Sq))

AvgSq =	zeros(Float64,length(q))
AvgRq = zeros(Float64,length(q))
alpha_q = zeros(Float64,length(q))

for j in 1:length(q)
	print("mS",q[j], ",mR", q[j], ",alpha", q[j], ",")
end

println()

for j in 1:length(q)
	AvgSq[j] = NaNMath.mean(filter(!ismissing,Sq[j]))
	AvgRq[j] = NaNMath.mean(filter(!ismissing,Rq[j]))
	alpha_q[j] = AvgSq[j] / (AvgRq[j] * log(lambda))
	print(AvgSq[j], ",", AvgRq[j], ",", alpha_q[j], ",")
end
println()

#s=sigma/100
#ww=W/100

#df = DataFrame(L=L,sig=s,w=ww,alpha_q1=alpha_q[1],alpha_q2=alpha_q[2],alpha_q3=alpha_q[3],alpha_q4=alpha_q[4],alpha_q5=alpha_q[5])

#CSV.open("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_10_21/results.csv")
#CSV.write("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_10_21/results.csv", df, header=false, delim = ",", append = true)
