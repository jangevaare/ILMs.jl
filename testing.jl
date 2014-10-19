include("population.jl")
include("events.jl")
include("inference.jl")

using Distributions, MCMC, PDMats

pop_db1 = pop_db_fun(625, MvNormal(eye(2).*2))

dist_mat1 = distance_mat_fun(pop_db1)

edb = event_db_fun(pop_db1)
random_infect(edb, 2, 1)
for t = 2:10
  infect_fun(dist_mat1, edb, convert(Float64,t), 1, 15)
  recover_fun(edb, convert(Float64,t), 0.1)
end

recovery_times1 = find_recovery_times(edb, true)

sa = susceptible_array_fun(edb, 10)
ia = infectious_array_fun(edb, 10)

sum(sa,1)
@time sir_loglikelihood=create_sir_loglikelihood(dist_mat1, sa, ia, recovery_times1)
@time sir_loglikelihood([1,5,0.1])

mymodel=model(sir_loglikelihood, init=[1,5,0.5])
@time chain1=run(mymodel, RWM(0.5), SerialMC(steps=1000, burnin=0))

describe(chain1)
