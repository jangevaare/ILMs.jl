using DataFrames, Distances, Distributions

function pop_db_fun(n, d)
  """
  Generate a population of size `n` distributed spatially according to
  a bivariate distribution, `d`
  """
  coordinates = rand(d, n)
  DataFrame(ind_id = 1:n, x = coordinates[:,1], y = coordinates[:,1])
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
