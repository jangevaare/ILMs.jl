function state_timeseries(event_db, cd="discrete")
	"""
	Count the number of individuals in each state at each event time
	"""
	susceptible = fill(0, length(event_db.event_times[event_db.event_times .< Inf]))
	infectious = fill(0, length(event_db.event_times[event_db.event_times .< Inf]))
	if(size(event_db.events)[2] == 4)
		recovered = fill(0, length(event_db.event_times[event_db.event_times .< Inf]))
		for t = 1:length(event_db.event_times[event_db.event_times .< Inf])
			susceptible[t] = sum(find_state(event_db, event_db.event_times[t], "S", cd))
			infectious[t] = sum(find_state(event_db, event_db.event_times[t], "I", cd))
			recovered[t] = sum(find_state(event_db, event_db.event_times[t], "R", cd))
		end
		return DataFrame(time = event_db.event_times[event_db.event_times .< Inf], S = susceptible, I = infectious, R = recovered)
	end
    if(size(event_db.events)[2] == 3)
        for t = 1:length(event_db.event_times[event_db.event_times .< Inf])
            susceptible[t] = sum(find_state(event_db, event_db.event_times[t], "S", cd))
            infectious[t] = sum(find_state(event_db, event_db.event_times[t], "I", cd))
        end
        return DataFrame(time = event_db.event_times[event_db.event_times .< Inf], S = susceptible, I = infectious)
    end
end
