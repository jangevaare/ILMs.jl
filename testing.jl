include("population.jl")
include("events.jl")
include("inference.jl")

using Distributions, PDMats, Mamba

pop_db1 = pop_db_fun(625, MvNormal(eye(2).*2))
dist_mat1 = distance_mat_fun(pop_db1)
@time dist_mat_ab1 = distance_mat_alphabeta_fun(dist_mat1, 1, 15)
edb = event_db_fun(pop_db1)
random_infect(edb, 2, 1)

for t = 2:10
  infect_fun(dist_mat_ab1, edb, convert(Float64,t))
  recover_fun(edb, convert(Float64,t), 0.1)
end

recovery_times1 = find_recovery_times(edb, true)
sa = susceptible_array_fun(edb, 10)
ia = infectious_array_fun(edb, 10)
sum(ia,1)

dist_mat_ab1
test1=infect_time_fun(dist_mat_ab1, find_infectious_fun(edb,2), find_susceptible_fun(edb,2))
minimum(test1)

rand(Exponential(test1[3]))

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


