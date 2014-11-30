type sdb
  susceptible_array::Array
  infectious_array::Array
  unique_event_times::Vector
  event_type::Vector
end

function state_array(event_db::edb)
  """
  Save the likelihood function from repetively determining state,
  by producing an array which contains this information for all 
  relevant time steps
  """
  sa=sdb(fill(false, (size(event_db.events)[1], length(unique(event_db.event_times[event_db.event_times<Inf])))), 
    fill(false, (size(event_db.events)[1], length(unique(event_db.event_times[event_db.event_times<Inf])))), 
    unique(event_db.event_times[event_db.event_times<Inf]), fill("I", length(unique(event_db.event_times[event_db.event_times<Inf]))))
  for i = 1:length(sa.unique_event_times)
    sa.susceptible_array[:,i] = find_state(event_db, sa.unique_event_times[i], "S")
    sa.infectious_array[:,i] = find_state(event_db, sa.unique_event_times[i], "I")
  end
  for i = 2:length(sa.unique_event_times)
    sa.susceptible_array[:,i] == sa.susceptible_array[:,i-1] && sa.event_type[i] = "R"
  end
  return sa
end

function find_recovery_times(event_db::edb, narm=true)
  """
  Determine recovery times for individuals
  """
  if narm == false
    return event_db.events[:,4]  - event_db.events[:,3]
  end
  if narm == true
    times = event_db.events[:,4]  - event_db.events[:,3]
    return times[isnan(times).==false]
  end
end


#=
function create_loglikelihood(distance_mat, event_db::edb)
  """
  Create a log likelihood function for continuous and discrete 
  SI and SIR models
  """
  if event_db.cd == "discrete" && size(event_db.events)[2]==3
    function loglikelihood(parameters)
      if parameters[1] <= 0 || parameters[2] <= 0
        return -Inf
      end
      dist_ab = dist_ab_mtx(distance_mat, parameters[1], parameters[2])
      event_loglikes = fill(-Inf, length(state_array.unique_event_times)


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
=#