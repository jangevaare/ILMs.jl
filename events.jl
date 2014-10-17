using DataFrames

function event_db_fun(pop_db)
"""
Generate an event database with a single initial infection
"""
  DataFrame(ind_id = [pop_db[:ind_id], sample(pop_db[:ind_id])], newstatus = [rep('s',pop_db[end,:ind_id]),'i'], time = [rep(0,pop_db[end,:ind_id]), 1.0])
end

function find_susceptible_fun(event_db, time)
"""
Find individuals which have not been infected prior to `time`
"""
  step1 = event_db[event_db[:time] .< time,:]
  susceptible = step1[step1[:newstatus] .== 's', 1]
  infected = step1[step1[:newstatus] .== 'i', 1]
if length(infected) > 0
  step2 = trues(length(susceptible))
  for i = 1:length(susceptible)
    step2[i] = any(infected .== susceptible[i]) && false
  end
  susceptible[step2]
  else
    susceptible
  end
end

function find_infectious_fun(event_db, time)
"""
Find individuals which have been infected, but not recovered prior to `time`
"""
  step1 = event_db[event_db[:time] .< time,:]
  infected = step1[step1[:newstatus] .== 'i', 1]
  recovered = step1[step1[:newstatus] .== 'r', 1]
if length(recovered) > 0
  step2 = trues(length(infected))
  for i = 1:length(infected)
    step2[i] = any(recovered .== infected[i]) && false
  end
  infected[step2]
  else
    infected
  end
end

function find_recovered_fun(event_db, time)
"""
Find individuals which have not been infected prior to `time`
"""
  step1=event_db[event_db[:time] .< time,:]
  step1[step1[:newstatus] .== 'r', 1]
end
