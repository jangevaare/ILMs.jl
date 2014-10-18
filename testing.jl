include("population.jl")
include("events.jl")
using Distributions

pop_db1 = pop_db_fun(100)
dist_mat1 = distance_mat_fun(pop_db1)

edb = event_db_fun(pop_db1)
random_infect(edb, 2, 1)
for t = 1:20
  infect_fun(dist_mat1, edb, convert(Float64,t), 1, 5)
  recover_fun(edb, convert(Float64,t), 0.25)
end

recovery_times1 = find_recovery_times(edb)
loglikelihood(Geometric(0.8), recovery_times1)
