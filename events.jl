using DataFrames

function event_db_fun(pop_db)
"""
Generate an event database, where all individuals are suscpetible at time 0
"""
  DataFrame(ind_id = pop_db[:,1], s = zeros(size(pop_db)[1]), i = nans(size(pop_db)[1]), r = nans(size(pop_db)[1]))
end

function random_infect(event_db)
"""
  Randomly infect 1 individual at time = 1.0. A random infection is required
  at the beginning of a simulation. Used for SI models
"""
  event_db[sample(event_db[:,1], 1), 3] = 1.0
  event_db = event_db[:, [1,2,3]]
end

function random_infectrecover(event_db, exp_gamma)
"""
  Randomly infect 1 individual at time = 1.0, and their recovery based on a mean
  recovery time of `exp_gamma`. At least one random infection is required
  at the beginning of a simulation
"""
  chosen_one=sample(event_db[:,1], 1)
  event_db[chosen_one, 3] = 1.0
  event_db[chosen_one, 4] = 1.0 + rand(Exponential(1/exp_gamma), 1)
end

function find_susceptible_fun(event_db, time)
"""
Find individuals which have not been infected prior to `time`
"""
  susceptible_index = falses(size(event_db)[1])
  for i = 1:length(susceptible_index)
    susceptible_index[i] = isnan(event_db[i,3]) || event_db[i,3] >= time
  end
  susceptible_index
end

function find_infectious_fun(event_db, time)
"""
Find individuals which have been infected prior to `time`,
but have not recovered by `time`
"""
  infectious_index = falses(size(event_db)[1])
  for i = 1:length(infectious_index)
    infectious_index[i] = event_db[i,3] < time && ((isnan(event_db[i,4])) || (time <= event_db[i,4]))
  end
  infectious_index
end

function find_recovered_fun(event_db, time)
"""
Find individuals which have recovered prior to or at `time`
"""
  recovered_index = falses(size(event_db)[1])
  for i = 1:length(recovered_index)
    recovered_index[i] = event_db[i,4] <= time
  end
  recovered_index
end

function find_now_susceptible_fun(event_db, time)
"""
Find individuals which have not been infected prior to or at `time`
"""
  susceptible_index = falses(size(event_db)[1])
  for i = 1:length(susceptible_index)
    susceptible_index[i] = isnan(event_db[i,3]) || event_db[i,3] > time
  end
  susceptible_index
end

function find_now_infectious_fun(event_db, time)
"""
Find individuals which have been infected prior to or at `time`, and have either
not recovered, or have yet to recover
"""
  infectious_index = falses(size(event_db)[1])
  for i = 1:length(infectious_index)
    infectious_index[i] = event_db[i,3] <= time && ((isnan(event_db[i,4])) || (time < event_db[i,4]))
  end
  infectious_index
end

function find_recovery_times(event_db, narm)
"""
Determine recovery times for individuals
"""
  if narm == false
    event_db[:,4]  - event_db[:,3]
  elseif narm == true
    times = event_db[:,4]  - event_db[:,3]
    times[isnan(times).==false]
  else
    error("narm (2nd argument) must be boolean")
  end
end

function distance_mat_alphabeta_fun(distance_mat, alpha, beta)
"""
Find the -alpha*(distance matrix^-beta)
"""
  alpha .* (distance_mat .^ -beta)
end

function infect_prob_fun(distance_mat_alphabeta, infectious, susceptible)
"""
Determine infection probabilities for all susceptible individuals
"""
  1 .- exp(-sum(distance_mat_alphabeta[susceptible, infectious], 2))
end

function infect_fun(distance_mat_alphabeta, event_db, time)
"""
Propagate infection through population according to `alpha` and `beta`
"""
  susceptible = find_susceptible_fun(event_db, time)
  infect_probs = zeros(length(susceptible))
  infect_probs[susceptible]=infect_prob_fun(distance_mat_alphabeta, find_infectious_fun(event_db, time), susceptible)
  infected = falses(length(susceptible))
  for i in 1:length(infect_probs)
    infected[i] = rand() < infect_probs[i]
  end
  event_db[infected, 3] = time
end

function infect_time_fun(distance_mat_alphabeta, infectious, susceptible)
"""
Generate infection times (exponentially distributed) based on current
infectious and susceptible. The minimum time will become the next
infected individual, remaining times will need to be recalculated
"""
  exponential_rate = fill(Inf, length(susceptible))
  exponential_rate[susceptible] = sum(distance_mat_alphabeta[susceptible, infectious], 2) .^ -1
  infect_times = fill(Inf, length(susceptible))
  for i = 1:length(susceptible)
    infect_times[i]= rand(Exponential(exponential_rate[i]))
  end
  infect_times
end

function continuous_infect_fun(event_db, distance_mat_alphabeta)
"""
Generate an additional infection the lastest state of the population
"""
  maxtime_infectious=maximum(event_db[:,3])
  susceptible = find_now_susceptible_fun(event_db, maxtime_infectious)
  infectious = find_now_infectious_fun(event_db, maxtime_infectious)
  infection_times = infect_time_fun(distance_mat_alphabeta, infectious, susceptible)
  which_min = infection_times .== minimum(infection_times)
  event_db[which_min,3] = maxtime_infectious + infection_times(which_min)
end

function recover_fun(event_db, time, gamma_inverse)
"""
Recover individuals following a geometric distribution with p = `gamma_inverse`
"""
  infectious = find_infectious_fun(event_db, time)
  recovered = falses(length(infectious))
  for i = 1:length(infectious)
    recovered[i] = infectious[i] && (rand() < gamma_inverse)
  end
  event_db[recovered, 4] = rep(time, sum(recovered))
end
