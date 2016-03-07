"""
Count the number of individuals in each state at each event time
"""
function state_timeseries(event_db::edb)
	susceptible = fill(0, length(event_db.event_times[event_db.event_times .< Inf]))
	infectious = fill(0, length(event_db.event_times[event_db.event_times .< Inf]))
	if size(event_db.events)[2] == 4
		recovered = fill(0, length(event_db.event_times[event_db.event_times .< Inf]))
		for t = 1:length(event_db.event_times[event_db.event_times .< Inf])
			susceptible[t] = sum(find_state(event_db, event_db.event_times[t], "S"))
			infectious[t] = sum(find_state(event_db, event_db.event_times[t], "I"))
			recovered[t] = sum(find_state(event_db, event_db.event_times[t], "R"))
		end
		return DataFrame(time = event_db.event_times[event_db.event_times .< Inf], S = susceptible, I = infectious, R = recovered)
	end
  if size(event_db.events)[2] == 3
    for t = 1:length(event_db.event_times[event_db.event_times .< Inf])
        susceptible[t] = sum(find_state(event_db, event_db.event_times[t], "S"))
        infectious[t] = sum(find_state(event_db, event_db.event_times[t], "I"))
    end
    return DataFrame(time = event_db.event_times[event_db.event_times .< Inf], S = susceptible, I = infectious)
  end
end


"""
Determine infection status for each individual at a specific time,
and create a dataframe which includes their coordinates for easy
spatial plotting of status through time
"""
function spatial_states(pop_db::DataFrame, event_db::edb, time::Real)
  i = find_state(event_db, time, "I")
  state = fill("S", length(i))
  state[i] = "I"
  if size(event_db.events)[2] == 4
      state[find_state(event_db, time, "R")] = "R"
  end
  return append!(DataFrame(x=fill(NaN, 3), y=fill(NaN, 3), state=["S", "I", "R"]), DataFrame(x=pop_db[:,2], y=pop_db[:,3], state=state))
end


"""
Make a spatial animation/slider in iJulia
"""
function state_animation(framewindow, event_db, pop_db)
  @manipulate for time = 1:framewindow:ceil(maximum(event_db.event_times[event_db.event_times .< Inf]))
      plot(spatial_states(pop_db, event_db, time), x="x", y="y", Geom.point, color="state", Scale.discrete_color_manual("black", "orange", "blue"))
  end
end
