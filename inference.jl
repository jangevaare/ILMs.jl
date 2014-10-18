include("events.jl")
using Distributions

function SIR_loglikelihood(distance_mat, recovery_times, event_db, alpha, beta, gamma_inverse)
  loglikelihood(Geometric(gamma_inverse), recovery_times) +
  daily_loglikes = zeros(convert(int64, event_db[end, :time])-1)
  for i = 1:length(daily_loglikes)
    daily_loglikes[i] = infect_prob_fun(distance_mat, infectious, susceptible, alpha, beta)
