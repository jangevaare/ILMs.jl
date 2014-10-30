include("population.jl")
include("events.jl")

#include("inference.jl")
#using Distributions, PDMats, Mamba
#using DataFrames, Distances, Distributions

pop_db1 = create_pop_db(50, MvNormal(eye(2).*2))
dist_mat1 = create_dist_mtx(pop_db1)
dist_mat_ab1 = dist_ab_mtx(dist_mat1, 1, 15)

edb = create_event_db(pop_db1, "SIR")

intial_infect(edb, "continuous", 5)

recovery_times1 = find_recovery_times(edb, true)
sa = susceptible_array_fun(edb, 10)
ia = infectious_array_fun(edb, 10)
sum(ia,1)

sir_loglikelihood=create_sir_loglikelihood(dist_mat1, sa, ia, recovery_times1)
@time sir_loglikelihood([1,5,0.5])

n = 10000
burnin = 5000
sim = Chains(n, 3, names = ["alpha", "beta", "inversegamma"])
theta = AMMVariate([1.0, 5.0, 0.2])
SigmaF = cholfact(eye(3))
@time for i in 1:n
  amm!(theta, SigmaF, sir_loglikelihood, adapt = (i <= burnin))
  sim[i,:,1] = [theta[1:2], exp(theta[3])]
end
Mamba.describe(sim)
plot(sim)
