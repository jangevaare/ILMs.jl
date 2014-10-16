using DataFrames

function event_db_fun(pop_db)
"""
Generate an event database with a single initial infection
"""
  DataFrame(ind_id = sample(pop_db[:ind_id]), newstatus = ['i'], time = 1.0)
end

function find_infectious_fun(event_db, time)
"""
Find individuals which have been infected prior to `time`
"""
  event_db[event_db[:time] .< time,][event_db[:newstatus] .== 'i', :ind_id]
end

