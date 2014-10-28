include("population.jl")
include("events.jl")
include("inference.jl")

using Distributions, PDMats, Mamba

pop_db1 = pop_db_fun(50, MvNormal(eye(2).*2))
dist_mat1 = distance_mat_fun(pop_db1)
dist_mat_ab1 = distance_mat_alphabeta_fun(dist_mat1, 1, 15)

edb = event_db_fun(pop_db1)
random_infect(edb)
edb

edb = event_db_fun(pop_db1)
random_infectrecover(edb, 4)
edb

continuous_infect_fun(edb, dist_mat_ab1)

edb

([0.1:0.1:0.5])[[true, true, false, false, true]]

for t = 1:20
  continuous_infect_fun(edb, dist_mat_ab1)
end

print(maximum(edb[:,3])) + minimum(infection_times1))

sum(find_now_susceptible_fun(edb, 1.0))
sum(find_now_infectious_fun(edb, 1.0))

infection_times1 = infect_time_fun(dist_mat_ab1, find_now_infectious_fun(edb, 1.0), find_now_susceptible_fun(edb, 1.0))
infection_times1[infection_times1 .== minimum(infection_times1)]




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


