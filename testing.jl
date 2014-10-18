include("population.jl")
include("events.jl")

pop_db1 = pop_db_fun(5000)
dist_mat1 = distance_mat_fun(pop_db1)
edb = event_db_fun(pop_db1)
for t = 1:20
  infect_fun(dist_mat1, edb, convert(Float64,t), 1, 5)
  recover_fun(edb, convert(Float64,t), 0.25)
end


edb[edb[:ind_id] .== 400,]


find_recovery_times(edb)
