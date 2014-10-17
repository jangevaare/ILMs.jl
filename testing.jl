include("population.jl")
include("events.jl")

pop_db1 = pop_db_fun(500)

dist_mat1 = distance_mat_fun(pop_db1)

edb = event_db_fun(pop_db1)

find_infectious_fun(edb, 0.)
