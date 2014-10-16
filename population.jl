using DataFrames
using Distances

function pop_db_fun(n)
  """
  Generate a population of size `n` distributed spatially according to
  a bivariate normal distribution
  """
  cbind(DataFrame(ind_id = 1:n), DataFrame(randn(n, 2)))
end

function distance_mat_fun(pop_db)
  """
  Construct a matrix containing the euclidean distance between every
  pair of individuals
  """
  for i = pop_db[:ind_id]
    for j = pop_db[:ind_id]
      euclidean(pop_db[i, 2:3], pop_db[j, 2:3])


pop_db1 = pop_db_fun(100)
sample(1:100)

      euclidean(Vector(pop_db1[1, 2:3]), Vector(pop_db1[2, 2:3]))
      pop_db1[1,1:3]

DataFrame(ind_id = sample(1:100), newstatus = ['i'])
