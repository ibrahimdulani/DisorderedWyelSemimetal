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

            nboot=100                #number of bootstrap samples
            nsample=100              #number of data points in each bootstrap sample

            #average on bootstrap samples
            
            Sq_b = zeros(nboot,length(q))
            Rq_b = zeros(nboot,length(q))
            alpha_q_b = zeros(nboot,length(q))

            for k in 1:nboot

                for j in 1:length(q)
                    Sq_b[k,j] = NaNMath.mean(sample(Sq[j],nsample))
                    Rq_b[k,j] = NaNMath.mean(sample(Rq[j],nsample))
                    alpha_q_b[k,j] = Sq_b[k,j]/(Rq_b[k,j]*log(lambda))
                end

            end

            #bootstrap error 

            alpha = zeros(length(q))
            err_alpha = zeros(length(q))

            println()

            for jj in 1:length(q)
                alpha[jj] = mean(alpha_q_b[:,jj])
                err_alpha[jj] = std(alpha_q_b[:,jj])
                print(alpha[jj], ",", err_alpha[jj], ",")
            end

            println()

            s=sigma/100
			ww=W/100

            df = DataFrame(L=L,sig=s,w=ww,alpha_q0=alpha[1], err0=err_alpha[1], alpha_q04=alpha[2], err04=err_alpha[2], alpha_q05=alpha[3], err05=err_alpha[3], 
            alpha_q06=alpha[4], err06=err_alpha[4], alpha_q08=alpha[5], err08=err_alpha[5], alpha_q15=alpha[6],err15=err_alpha[6])

            CSV.open("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/results.csv")
			CSV.write("/Users/ibrahim/Documents/Rutgers/DisorderedWeylSemimetals/from_cluster/mw936/data_12_21/results.csv", df, header=false, delim = ",", append = true)

        end

    end

end