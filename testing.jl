include("population.jl")
include("events.jl")

pop_db1 = pop_db_fun(500)

dist_mat1 = distance_mat_fun(pop_db1)

edb = event_db_fun(pop_db1)

test1 = find_infectious_fun(edb, 2)
test2 = find_susceptible_fun(edb, 2)
test3 = find_recovered_fun(edb, 2)

infectious1 = find_infectious_fun(edb, 2)
susceptible1 = find_susceptible_fun(edb, 2)
infect_prob1=infect_prob_fun(dist_mat1, infectious1, susceptible1, 1, 8)

infect_fun(dist_mat1, edb, 2., 1, 5)
edb

