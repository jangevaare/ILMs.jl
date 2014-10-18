include("population.jl")
include("events.jl")
include("inference.jl")
using Distributions

pop_db1 = pop_db_fun(100)
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

@time SIR_loglikelihood(dist_mat1, sa, ia, recovery_times1, 1, 5, 0.1)
