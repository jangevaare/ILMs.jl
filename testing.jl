include("population.jl")
include("events.jl")
include("inference.jl")

using Distributions, PDMats

pop_db1 = pop_db_fun(1000)
dist_mat1 = distance_mat_fun(pop_db1)

edb = event_db_fun(pop_db1)
random_infect(edb, 2, 1)
for t = 2:20
  infect_fun(dist_mat1, edb, convert(Float64,t), 1, 5)
  recover_fun(edb, convert(Float64,t), 0.1)
end

recovery_times1 = find_recovery_times(edb, true)

sa = susceptible_array_fun(edb, 20)
ia = infectious_array_fun(edb, 20)

rand(MvNormal(PDMat([[1 0], [0 1]])), 10)

@time sir_loglikelihood=create_sir_loglikelihood(dist_mat1, sa, ia, recovery_times1)

@time sir_loglikelihood(1,5,0.24)
