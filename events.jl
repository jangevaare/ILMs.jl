using DataFrames

event_db_fun(pop_db)
"""
Generate an event database with a single initial infection
"""
  DataFrame(ind_id = sample(pop_db[:ind_id]), newstatus = ['i'])
end

