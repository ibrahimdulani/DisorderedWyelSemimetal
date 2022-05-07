using Printf
using Statistics
using CSV
using DataFrames
using NaNMath        #package for ignoring nan values
using StatsBase

LL = [100, 500, 1000, 5000, 10000]
Ws = [1, 10, 15, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50]
sigmas = [1, 5, 10, 15, 20, 25, 30, 33, 35, 37, 40, 45, 46, 47, 48, 49, 50]
q = [0.0, 0.4, 0.5, 0.6, 0.8, 1.5]    #this should be the same with Sq_Rq.jl
lambda = 1/25     #this also should be the same with Sq_Rq.jl

for sigma in sigmas 

	for L in LL

		for W in Ws

			folder_name = @sprintf "/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/sigma%03d/L%d/W%03d/" sigma L W

			file_list = readdir(folder_name, join=true)

			file_list = filter!(x->stat(x).size>10,file_list)

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
#				filter!(ismissing,Sq[j])           #filters out missing values from Sq[j]
#				filter!(ismissing,Rq[j])           #filters out missing values from Rq[j]
				AvgSq[j] = NaNMath.mean(Sq[j])
				AvgRq[j] = NaNMath.mean(Rq[j])
				alpha_q[j] = AvgSq[j] / (AvgRq[j] * log(lambda))
				print(AvgSq[j], ",", AvgRq[j], ",", alpha_q[j], ",")
			end
			println()

			s=sigma/100
			ww=W/100

			df = DataFrame(L=L,sig=s,w=ww,alpha_q0=alpha_q[1], alpha_q04=alpha_q[2], alpha_q05=alpha_q[3],alpha_q06=alpha_q[4], alpha_q08=alpha_q[5],alpha_q15=alpha_q[6])

			CSV.open("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/results.csv")
			CSV.write("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/results.csv", df, header=false, delim = ",", append = true)

			print(" ", ww)
		end
	end
end