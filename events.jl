using DataFrames

function event_db_fun(pop_db)
"""
Generate an event database, where all individuals are suscpetible at time 0
"""
  DataFrame(ind_id = pop_db[:ind_id], s = zeros(size(pop_db)[1]), i = nans(size(pop_db)[1]), r = nans(size(pop_db)[1]))
end

function random_infect(event_db, n, time)
"""
Randomly infect `n` individuals at `time`. At least one random infection is required
at the beginning of a simulation
"""
  event_db[sample(event_db[:ind_id], n), 3] = time
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
Find individuals which have been infected prior to `time`, but have not recovered
"""
  infectious_index = falses(size(event_db)[1])
  for i = 1:length(infectious_index)
    infectious_index[i] = event_db[i,3] < time && ((isnan(event_db[i,4])) || (time <= event_db[i,4]))
  end
  infectious_index
end

function find_recovered_fun(event_db, time)
"""
Find individuals which have recovered prior to `time`
"""
  recovered_index = falses(size(event_db)[1])
  for i = 1:length(recovered_index)
    recovered_index[i] = event_db[i,4] < time
  end
  recovered_index
end

function find_recovery_times(event_db)
"""
Determine recovery times for individuals
"""
event_db[:,4]  - event_db[:,3]
end

function infect_prob_fun(distance_mat, infectious, susceptible, alpha, beta)
"""
Determine infection probabilities for all susceptible individuals
"""
  1 .- exp(-alpha .* sum(distance_mat[susceptible, infectious].^-beta, 2))
end

function infect_fun(distance_mat, event_db, time, alpha, beta)
"""
Propagate infection through population according to `alpha` and `beta`
"""
  susceptible = find_susceptible_fun(event_db, time)
  infect_probs=infect_prob_fun(distance_mat, find_infectious_fun(event_db, time), susceptible, alpha, beta)
  infected = falses(length(susceptible))
  for i in 1:length(susceptible)
    infected[i] = rand() < infect_probs[i,1]
  end
  event_db[infected, 3] = time
end

function recover_fun(event_db, time, gamma_inverse)
"""
Recover individuals following a geometric distribution with p = `gamma_inverse`
"""
  infectious = find_infectious_fun(event_db, time)
  recovered = falses(length(infectious))
  for i = 1:length(infectious)
    recovered[i] = infectious[i] && rand() < gamma_inverse
  end
  event_db[recovered, 4] = time
end
