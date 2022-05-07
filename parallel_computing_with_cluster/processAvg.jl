using Printf
using Statistics
using CSV
using DataFrames

LL = [100, 500, 1000, 5000, 10000]
Ws = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50]
sigmas = [5, 10, 20, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50]
q = [0.2, 0.5, 0.8, 1.0, 1.2]    #this should be the same with Sq_Rq.jl
lambda = 1/25     #this also should be the same with Sq_Rq.jl

for L in LL 

	for sigma in sigmas

		for W in Ws

			folder_name = @sprintf "/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data/L%d/sigma%03d/W%03d/" L sigma W

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
				AvgSq[j] = mean(Sq[j])
				AvgRq[j] = mean(Rq[j])
				alpha_q[j] = AvgSq[j] / (AvgRq[j] * log(lambda))
				print(AvgSq[j], ",", AvgRq[j], ",", alpha_q[j], ",")
			end
			println()

			s=sigma/100
			ww=W/100

			df = DataFrame(L=L,sig=s,w=ww,alpha_q1=alpha_q[1],alpha_q2=alpha_q[2],alpha_q3=alpha_q[3],alpha_q4=alpha_q[4],alpha_q5=alpha_q[5])

			CSV.open("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data/results.csv")
			CSV.write("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data/results.csv", df, header=false, delim = ",", append = true)

		end
	end
end
