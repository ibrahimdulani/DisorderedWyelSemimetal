using Printf
using Statistics
#using CSV
#using DataFrames

L = 100
W = 10
sigma = 10
q = [0.1, 0.4]    #this should be the same with Sq_Rq.jl
lambda = 1/25     #this also should be the same with Sq_Rq.jl

folder_name = @sprintf "/scratch/mw936/data/L%d/sigma%03d/W%03d/" L sigma W

file_list = readdir(folder_name, join=true)

print(size(file_list))
