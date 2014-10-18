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
    infectious_array[:,i] = find_infectious_fun(event_db, i+1)
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
    susceptible_array[:,i] = find_susceptible_fun(event_db, i+1)
  end
  susceptible_array
end

function create_sir_loglikelihood(distance_mat, susceptible_array, infectious_array, recovery_times)
"""
Create the SIR loglikelihood function, which accepts only 3 parameters,
alpha, beta, and gamma_inverse
"""
  function sir_loglikelihood(alpha, beta, gamma)
  """
  Compute the log likelihood of discrete time SIR models
  """
    daily_loglikes = zeros(sum(sum(susceptible_array,1) .> 0)-1)
    for i = 1:length(daily_loglikes)
      infect_probs=nans(size(susceptible_array)[1])
      infect_probs[susceptible_array[:,i]]=infect_prob_fun(distance_mat, infectious_array[:,i], susceptible_array[:,i], alpha, beta)
      daily_loglikes[i]=sum(log([(infect_probs[susceptible_array[:,i] & infectious_array[:,i+1]]), 1.-(infect_probs[susceptible_array[:,i] & susceptible_array[:,i+1]])]))
    end
    loglikelihood(Geometric(1/gamma), recovery_times) + sum(daily_loglikes)
  end
sir_loglikelihood
end

function mcmc_array_fun(iterations, parameters)
"""
Create a simple empty array for use in MCMC,
user must then insert initial parameters
"""
  zeros(iterations, parameters+1)
end

function mh_mcmc_fun(target, mcmc_array, prop_cov)
"""
Perform simple MH according to a log target and a previously
specified mcmc_array
"""
  i=1
  mcmc_array[i,end] = target(convert(Tuple, mcmc_array[i,1:(end-1)]))
  transitions = rand(MvNormal(prop_cov), size(mcmc_array)[1]-1)
  for i = 2:(size(mcmc_array)[1])

    target(convert(Tuple, mcmc_array[i-1,1:(end-1)]+transitions[i]))

