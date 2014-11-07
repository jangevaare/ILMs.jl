function state_timeseries(event_db, cd="discrete")
	"""
	Count the number of individuals in each state at each event time
	"""
	susceptible = fill(0, length(event_db.event_times))
	infectious = fill(0, length(event_db.event_times))
	if(size(event_db.events)[2] == 4)
		recovered = fill(0, length(event_db.event_times))
		for t = 1:length(event_db.event_times)
			susceptible[t] = sum(find_state(event_db, event_db.event_times[t], "S", cd))
			infectious[t] = sum(find_state(event_db, event_db.event_times[t], "I", cd))
			recovered[t] = sum(find_state(event_db, event_db.event_times[t], "R", cd))
		end
		return DataFrame(time = event_db.event_times, S = susceptible, I = infectious, R = recovered)
	end
    if(size(event_db.events)[2] == 3)
        recovered = fill(0, length(event_db.event_times))
        for t = 1:length(event_db.event_times)
            susceptible[t] = sum(find_state(event_db, event_db.event_times[t], "S", cd))
            infectious[t] = sum(find_state(event_db, event_db.event_times[t], "I", cd))
        end
        return DataFrame(time = event_db.event_times, S = susceptible, I = infectious)
    end
end
