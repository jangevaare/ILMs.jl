using DataFrames
using Distances

function pop_db_fun(n)
  """
  Generate a population of size `n` distributed spatially according to
  a bivariate normal distribution, more distribution settings may be
  defined in the future
  """
  cbind(DataFrame(ind_id = 1:n), DataFrame(randn(n, 2)))
end

function distance_mat_fun(pop_db)
  """
  Construct a matrix containing the euclidean distance between every
  pair of individuals contained in `pop_db`
  """
  distance_mat = zeros(pop_db[end, :ind_id],pop_db[end, :ind_id])
  for i = 1:pop_db[end,:ind_id]
    for j = (i+1):pop_db[end, :ind_id]
      distance_mat[j,i] = euclidean([pop_db[i, 2],pop_db[i, 3]], [pop_db[j, 2],pop_db[j, 3]])
    end
  end
  distance_mat' + distance_mat
end
