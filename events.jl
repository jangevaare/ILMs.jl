using DataFrames, Distributions

function create_event_db(pop_db, ilm="SIR")
  """
  Generate an event database, where all individuals are suscpetible at time
  0. Specify whether for an SI model or an SIR model.
  """
  if ilm == "SI"
    return DataFrame(ind_id = pop_db[:,1], s = fill(0, size(pop_db)[1]), i = fill(NaN, size(pop_db)[1]))
  end
  if ilm == "SIR"
    return DataFrame(ind_id = pop_db[:,1], s = fill(0, size(pop_db)[1]), i = fill(NaN, size(pop_db)[1]), r = fill(NaN, size(pop_db)[1]))
  end
end

function intial_infect(event_db, cd="discrete", gamma=0)
  """
  Randomly infect 1 individual at time = 1.0, A random infection is required
  at the beginning of all simulations. Specify whether ILM is continuous
  ("c") or discrete ("d"), which selects between a exponential or geometric
  recovery time functions. If `gamma` (the mean recovery time) is > 0, then
  a recovery for this infection will also be probabilistically generated.
  """
  if gamma > 0 && size(event_db)[2] == 4
    if cd == "continuous"
      recovery_dist = Exponential(gamma)
    end
    if cd == "discrete"
      recovery_dist = Geometric(1/gamma)
    end
    chosen_one=sample(event_db[:,1], 1)
    event_db[chosen_one, 3] = 1.0
    event_db[chosen_one, 4] = 1.0 + rand(recovery_dist, 1)
  else
    event_db[sample(event_db[:,1], 1), 3] = 1.0
  end
  return event_db
end

function find_state(event_db, time, state, cd="discrete")
  """
  Find the individuals falling into `state` (S, I, or R), at `time` from
  a continuous or discrete ILM
  """
  state_index = fill(false, size(event_db)[1])
  if state == "S"
    if cd == "discrete"
      for i = 1:length(state_index)
        state_index[i] = isnan(event_db[i,3]) || event_db[i,3] >= time
      end
    end
    if cd == "continuous"
      for i = 1:length(state_index)
        state_index[i] = isnan(event_db[i,3]) || event_db[i,3] > time
      end
    end
  end
  if state == "I"
    if cd == "discrete"
      for i = 1:length(state_index)
        state_index[i] = event_db[i,3] < time && ((isnan(event_db[i,4])) || (time <= event_db[i,4]))
      end
    end
    if cd == "continuous"
      for i = 1:length(state_index)
        state_index[i] = event_db[i,3] <= time && ((isnan(event_db[i,4])) || (time < event_db[i,4]))
      end
    end
  end
  if state == "R"
    if cd == "discrete"
      for i = 1:length(state_index)
        state_index[i] = event_db[i,4] < time
      end
    end
    if cd == "continuous"
      for i = 1:length(state_index)
        state_index[i] = event_db[i,4] <= time
      end
    end
  end
  return state_index
end

function find_recovery_times(event_db, narm=true)
  """
  Determine recovery times for individuals
  """
  if narm == false
    return event_db[:,4]  - event_db[:,3]
  elseif narm == true
    times = event_db[:,4]  - event_db[:,3]
    return times[isnan(times).==false]
  else
    error("narm (2nd argument) must be boolean")
  end
end

function dist_ab_mtx(distance_mat, alpha, beta)
  """
  Find the -alpha*(distance matrix^-beta)
  """
  return alpha .* (distance_mat .^ -beta)
end

function infection_probabilities(distance_mat_alphabeta, infectious, susceptible)
  """
  Determine infection probabilities for all susceptible individuals
  """
  1 .- exp(-sum(distance_mat_alphabeta[susceptible, infectious], 2))
end

function infection_times(distance_mat_alphabeta, infectious, susceptible)
  """
  Generate infection times (exponentially distributed) based on current
  infectious and susceptible. The minimum time will become the next
  infected individual, remaining times will need to be recalculated
  """
  exponential_rate = fill(Inf, length(susceptible))
  exponential_rate[susceptible] = sum(distance_mat_alphabeta[susceptible, infectious], 2)
  infect_times = fill(Inf, length(susceptible))
  for i = 1:length(susceptible)
    infect_times[i]= rand(Exponential(exponential_rate[i]))
  end
  infect_times
end

function infect_recover(distance_mat_alphabeta, event_db, cd="discrete", time=maximum(event_db[:,3]), gamma=0)
  """
  Propagate infection through a population, default is for a discrete model.
  Use of an SIR or SI model inferred from dimensions of `event_db` and
  specifcation of a `gamma` > 0. Default time of previous time step is based
  on the maximum infection time - which may not work well in discrete SIR
  models especially
  """
  susceptible = find_state(event_db, time, "S", cd)
  infectious = find_state(event_db, time, "I", cd)
  if cd == "discrete"
    infect_probs = fill(0, length(susceptible))
    infect_probs[susceptible]=infection_probabilities(distance_mat_alphabeta, infectious, susceptible)
    infected = fill(false, length(susceptible))
    for i in 1:length(infect_probs)
      infected[i] = rand() < infect_probs[i]
    end
    event_db[infected, 3] = time + 1
    if gamma > 0 && size(event_db)[2] == 4
      event_db[infected, 4] = (time + 2) .+ rand(Geometric(1/gamma), sum(infected))
    end
  end
  if cd == "continuous"
    infect_times = infection_times(distance_mat_alphabeta, infectious, susceptible)
    which_min = infect_times .== minimum(infect_times)
    event_db[which_min,3] = time + infect_times[which_min]
    if gamma > 0 && size(event_db)[2] == 4
      event_db[which_min,4] = event_db[which_min,3] + rand(Exponential(gamma),1)
    end
  end
  return event_db
end
