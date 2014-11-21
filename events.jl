type edb
  events::DataFrame
  event_times::Vector
  cd::ASCIIString
end

function create_event_db(pop_db, cd="discrete", gamma=Inf)
  """
  Generate an event database, where all individuals are suscpetible at time
  0. Specify whether for an SI model or an SIR model.
  """
  if gamma == Inf
    return edb(DataFrame(ind_id = pop_db[:,1], s = fill(0, size(pop_db)[1]), i = fill(NaN, size(pop_db)[1])), fill(Inf, size(pop_db)[1]), cd)
  end
  if 0 < gamma < Inf
    return edb(DataFrame(ind_id = pop_db[:,1], s = fill(0, size(pop_db)[1]), i = fill(NaN, size(pop_db)[1]), r = fill(NaN, size(pop_db)[1])), fill(Inf, 2*(size(pop_db)[1])), cd)
  end
end

function event_time_update!(eventtime, event_db)
  """
  Update the ordered sequence of event times based on a new event at
  `time`
  """
  maxevents=length(event_db.event_times)
  for i = 1:length(eventtime)
    event_db.event_times=[event_db.event_times[(1:maxevents)[eventtime[i] .>= event_db.event_times]], eventtime[i], event_db.event_times[(1:(maxevents-1))[eventtime[i] .< (event_db.event_times[1:(maxevents-1)])]]]
  end
end

function initial_infect!(event_db, gamma=Inf)
  """
  Randomly infect 1 individual at time = 1.0, A random infection is required
  at the beginning of all simulations. Specify whether ILM is continuous
  ("c") or discrete ("d"), which selects between a exponential or geometric
  recovery time functions. If `gamma` (the mean recovery time) is > 0, then
  a recovery for this infection will also be probabilistically generated.
  """
  chosen_one=sample(event_db.events[:,1], 1)
  event_db.events[chosen_one, 3] = 1.0
  event_time_update(1.0, event_db)
  if 0 < gamma < Inf && size(event_db.events)[2] == 4
    if event_db.cd == "continuous"
      recovery_dist = Exponential(gamma)
    end
    if event_db.cd == "discrete"
      recovery_dist = Geometric(1/gamma)
    end
    recovery_time = 1.0 + rand(recovery_dist, 1)
    event_db.events[chosen_one, 4] = recovery_time
    event_time_update(recovery_time, event_db)
  end
end

function find_state(event_db, time, state)
  """
  Find the individuals falling into `state` (S, I, or R), at `time` from
  a continuous or discrete ILM
  """
  state_index = fill(false, size(event_db.events)[1])
  if state == "S"
    if event_db.cd == "discrete"
      for i = 1:length(state_index)
        state_index[i] = isnan(event_db.events[i,3]) || event_db.events[i,3] >= time
      end
    end
    if event_db.cd == "continuous"
      for i = 1:length(state_index)
        state_index[i] = isnan(event_db.events[i,3]) || event_db.events[i,3] > time
      end
    end
  end
  if state == "I"
    if size(event_db.events)[2] == 4
      if event_db.cd == "discrete"
        for i = 1:length(state_index)
          state_index[i] = event_db.events[i,3] < time && ((isnan(event_db.events[i,4])) || (time <= event_db.events[i,4]))
        end
      end
      if event_db.cd == "continuous"
        for i = 1:length(state_index)
          state_index[i] = event_db.events[i,3] <= time && ((isnan(event_db.events[i,4])) || (time < event_db.events[i,4]))
        end
      end
    end
    if size(event_db.events)[2] == 3
      if event_db.cd == "discrete"
        for i = 1:length(state_index)
          state_index[i] = event_db.events[i,3] < time
        end
      end
      if event_db.cd == "continuous"
        for i = 1:length(state_index)
          state_index[i] = event_db.events[i,3] <= time
        end
      end
    end
  end
  if state == "R"
    if event_db.cd == "discrete"
      for i = 1:length(state_index)
        state_index[i] = event_db.events[i,4] < time
      end
    end
    if event_db.cd == "continuous"
      for i = 1:length(state_index)
        state_index[i] = event_db.events[i,4] <= time
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
    return event_db.events[:,4]  - event_db.events[:,3]
  elseif narm == true
    times = event_db.events[:,4]  - event_db.events[:,3]
    return times[isnan(times).==false]
  else
    error("narm (2nd argument) must be boolean")
  end
end

function dist_ab_mtx(distance_mat, alpha, beta)
  """
  Find the alpha*(distance matrix^-beta)
  """
  return alpha .* (distance_mat .^ -beta)
end

function infection_probabilities(distance_mat_alphabeta, infectious, susceptible)
  """
  Determine infection probabilities for all susceptible individuals
  """
  return 1 .- exp(-sum(distance_mat_alphabeta[susceptible, infectious], 2))
end

function infection_times(distance_mat_alphabeta, infectious, susceptible)
  """
  Generate infection times (exponentially distributed) based on current
  infectious and susceptible. The minimum time will become the next
  infected individual, remaining times will need to be recalculated
  """
  exponential_rate = fill(Inf, length(susceptible))
  exponential_rate[susceptible] = sum(distance_mat_alphabeta[susceptible, infectious], 2).^-1.0
  infect_times = fill(Inf, length(susceptible))
  for i = 1:length(susceptible)
    if 0 < exponential_rate[i] < Inf 
      infect_times[i]= rand(Exponential(exponential_rate[i]))
    end
  end
  return infect_times
end

function infect_recover!(distance_mat_alphabeta, event_db, time=1.0, gamma=Inf)
  """
  Propagate infection through a population, default is for a discrete model.
  Use of an SIR or SI model inferred from dimensions of `event_db` and
  specifcation of a `gamma` > 0. Default time of previous time step is based
  on the maximum infection time - which may not work well in discrete SIR
  models especially
  """
  susceptible = find_state(event_db, time, "S")
  infectious = find_state(event_db, time, "I")
  if event_db.cd == "discrete"
    infect_probs = fill(0, length(susceptible))
    infect_probs[susceptible]=infection_probabilities(distance_mat_alphabeta, infectious, susceptible)
    infected = fill(false, length(susceptible))
    for i in 1:length(infect_probs)
      infected[i] = rand() < infect_probs[i]
    end
    event_db[infected, 3] = time + 1
    if 0 < gamma < Inf && size(event_db.events)[2] == 4
      recovery_times = rand(Geometric(1/gamma), sum(infected))
      event_db.events[infected, 4] = (time + 2) .+ recovery_times
      for i = 1:sum(infected)
        event_time_update!(time+1, event_db)
        event_time_update!(time+2+recovery_times[i], event_db)
      end
    end
  end
  if event_db.cd == "continuous"
    if gamma == Inf && size(event_db.events)[2] == 3
      infect_times = infection_times(distance_mat_alphabeta, infectious, susceptible)
      which_min = infect_times .== minimum(infect_times)
      event_db.events[which_min,3] = time + infect_times[which_min]
      event_time_update!(time + infect_times[which_min], event_db)
    end
    if 0 < gamma < Inf && size(event_db.events)[2] == 4
      infect_times = infection_times(distance_mat_alphabeta, infectious, susceptible)
      which_min = infect_times .== minimum(infect_times)
      while (time .< event_db.event_times) != ((time + infect_times[which_min][1]) .<= event_db.event_times)
        time = event_db.event_times[time .< event_db.event_times][1]
        susceptible = find_state(event_db, time, "S")
        infectious = find_state(event_db, time, "I")
        infect_times = infection_times(distance_mat_alphabeta, infectious, susceptible)
        if minimum(infect_times) == Inf
          break
        end
        which_min = infect_times .== minimum(infect_times)
      end
      if minimum(infect_times) != Inf
        event_db.events[which_min,3] = (time + infect_times[which_min])
        recovery_time = rand(Exponential(gamma),sum(which_min))
        event_db.events[which_min,4] = event_db.events[which_min,3] + recovery_time
        event_time_update!(time + infect_times[which_min], event_db)
        event_time_update!(time + infect_times[which_min] + recovery_time, event_db)
      end
    end
  end
end

function infect_recover_loop(pop_db, cd="discrete", ilm="SI", alpha=1, beta=1, gamma=Inf)
  """
  Function to generate initial infection then loop infect_recover function as appropriate 
  for continuous and discrete SI and SIR models
  """
  distance_mat_alphabeta = dist_ab_mtx(create_dist_mtx(pop_db), alpha, beta)
  event_db = create_event_db(pop_db, cd, gamma)
  initial_infect!(event_db, gamma)
  time=1.0
  if event_db.cd == "discrete"
    if gamma == Inf && size(event_db.events)[2] == 3
      while sum(find_state(event_db, time, "S")) > 0
        infect_recover!(distance_mat_alphabeta, event_db, time, Inf)
        time = time + 1.0
      end
    end
    if 0 < gamma < Inf && size(event_db.events)[2] == 4
      while sum(find_state(event_db, time, "I")) > 0
        infect_recover!(distance_mat_alphabeta, event_db, time, gamma)
        time = time + 1.0
      end
    end
  end
  if event_db.cd == "continuous"
    if gamma == Inf && size(event_db.events)[2] == 3
      while sum(find_state(event_db, time, "S")) > 0 && time < Inf
        infect_recover!(distance_mat_alphabeta, event_db, time, Inf)
        time = maximum(event_db.events[:,3])
      end
    end
    if 0 < gamma < Inf && size(event_db.events)[2] == 4
      while sum(find_state(event_db, time, "I")) > 0 && time < Inf
        infect_recover!(distance_mat_alphabeta, event_db, time, gamma)
        if time == maximum(event_db.events[:,3])
          break
        end
        time = maximum(event_db.events[:,3])
      end
    end
  end
  return event_db
end