"""
Generate a population of size `n` distributed spatially according to
a bivariate distribution, `d`
"""
function create_pop_db_bivariate(n::Int64, d::ContinuousMultivariateDistribution)
  coordinates = rand(d, n)'
  DataFrame(ind_id = 1:n, x = coordinates[:,1], y = coordinates[:,2])
end


"""
Generate a population of size `n` distributed spatially according to
a univariate distributions, `d1` and `d2`
"""
function create_pop_db_univariates(n::Int64,
                                   d1::ContinuousUnivariateDistribution,
                                   d2::ContinuousUnivariateDistribution)
  DataFrame(ind_id = 1:n, x = rand(d1, n), y = rand(d2, n))
end


"""
Construct a matrix containing the euclidean distance between every
pair of individuals contained in `pop_db`
"""
function create_dist_mtx(pop_db::DataFrame)
  distance_mat = fill(0.0, (pop_db[end, 1], pop_db[end, 1]))
  for i = 1:(size(distance_mat)[1]-1)
    for j = (i+1):(size(distance_mat)[1])
      distance_mat[i,j] = euclidean([pop_db[i, 2],pop_db[i, 3]], [pop_db[j, 2],pop_db[j, 3]])
    end
  end
  distance_mat' + distance_mat
end
