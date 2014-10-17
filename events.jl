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
  step1=event_db[event_db[:time] .< time,:]
  susceptible=step1[step1[:newstatus] .== 's', 1]
  step2 = step1[step1[:newstatus] .== 'i', 1]
if length(step2) > 0
  step3 = trues(length(step1))
  for i = 1:length(step1)
    step3[i] = any(step2 .== step1[i]) && false
  end
  susceptible[step3]
  else
    susceptible
  end
end

function find_infectious_fun(event_db, time)
"""
Find individuals which have been infected, but not recovered prior to `time`
"""
  step1=event_db[event_db[:time] .< time,:]
  infected=step1[step1[:newstatus] .== 'i', 1]
  step2 = step1[step1[:newstatus] .== 'r', 1]
if length(step2) > 0
  step3 = trues(length(step1))
  for i = 1:length(step1)
    step3[i] = any(step2 .== step1[i]) && false
  end
  infected[step3]
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
