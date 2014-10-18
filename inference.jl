include("events.jl")
using Distributions

function infectious_array_fun(event_db, obs_length)
"""
Save the liklihood function from repetively determining infectious status,
by producing an array which contains this information for all relevant
time steps
"""
  infectious_array=falses(size(event_db)[1], obs_length)
  for i = 1:obs_length
    infectious_array[:,i] = find_infectious_fun(event_db, i)
  end
  infectious_array
end

function susceptible_array_fun(event_db, obs_length)
"""
Save the liklihood function from repetively determining susceptible status,
by producing an array which contains this information for all relevant
time steps
"""
  susceptible_array=falses(size(event_db)[1], obs_length)
  for i = 1:obs_length
    susceptible_array[:,i] = find_susceptible_fun(event_db, i)
  end
  susceptible_array
end

function SIR_loglikelihood(distance_mat, recovery_times, event_db, alpha, beta, gamma_inverse, obs_length)
  loglikelihood(Geometric(gamma_inverse), recovery_times) +
  daily_loglikes = zeros(obs_length-1)
  for i = 1:length(daily_loglikes)
    infectious = find_infectious_fun
    daily_loglikes[i] = infect_prob_fun(distance_mat, infectious, susceptible, alpha, beta)
