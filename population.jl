using DataFrames

function pop_db(n)
  """
  Generate a population of size `n` distributed spatially according to
  a bivariate normal distribution
  """
  cbind(DataFrame(ind_id = 1:n), DataFrame(randn(n, 2)))
end
