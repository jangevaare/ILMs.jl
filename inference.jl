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
  function sir_loglikelihood(parameters)
  """
  Compute the log likelihood of discrete time SIR models (alpha, beta, gamma)
  """
    distance_mat_alphabeta = distance_mat_alphabeta_fun(distance_mat, parameters[1], parameters[2])
    if parameters[1] <= 0 || parameters[2] <= 0 || parameters[3] >= 1 || parameters[3] <= 0
      -Inf
    else
      daily_loglikes = zeros(sum(sum(susceptible_array,1) .> 0)-1)
      for i = 1:length(daily_loglikes)
        infect_probs=nans(size(susceptible_array)[1])
        infect_probs[susceptible_array[:,i]]=infect_prob_fun(distance_mat_alphabeta, infectious_array[:,i], susceptible_array[:,i])
        daily_loglikes[i]=sum(log([(infect_probs[susceptible_array[:,i] & infectious_array[:,i+1]]), 1.-(infect_probs[susceptible_array[:,i] & susceptible_array[:,i+1]])]))
      end
      loglikelihood(Geometric(parameters[3]), recovery_times) + logpdf(Uniform(0,1), parameters[3]) + logpdf(Gamma(1,1), parameters[1]) + logpdf(Gamma(5,1), parameters[2]) + sum(daily_loglikes)
    end
  end
sir_loglikelihood
end
