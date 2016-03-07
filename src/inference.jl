type sdb
  susceptible::Array{Bool,2}
  infectious::Array{Bool,2}
  recovered::Array{Bool,2}
  new_infection::Array{Bool,2}
  new_recovery::Array{Bool,2}
  event_times::Array{Float64,1}
  event_type::Array{ASCIIString,1}
end


"""
Save the likelihood function from repetively determining state,
by producing an array which contains this information for all
relevant time steps
"""
function state_array(event_db::edb)
  sa=sdb(fill(false, (size(event_db.events)[1], length(unique(event_db.event_times[event_db.event_times.<Inf])))), fill(false, (size(event_db.events)[1], length(unique(event_db.event_times[event_db.event_times.<Inf])))), fill(false, (size(event_db.events)[1], length(unique(event_db.event_times[event_db.event_times.<Inf])))), fill(false, (size(event_db.events)[1], length(unique(event_db.event_times[event_db.event_times.<Inf])))), fill(false, (size(event_db.events)[1], length(unique(event_db.event_times[event_db.event_times.<Inf])))), unique(event_db.event_times[event_db.event_times.<Inf]), fill("I", length(unique(event_db.event_times[event_db.event_times.<Inf]))))
  for i = 1:length(sa.event_times)
    sa.susceptible[:,i] = find_state(event_db, sa.event_times[i], "S")
    sa.infectious[:,i] = find_state(event_db, sa.event_times[i], "I")
    if size(event_db.events)[2]==4
      sa.recovered[:,i] = find_state(event_db, sa.event_times[i], "R")
    end
  end
  for i = 2:length(sa.event_times)
    if sa.susceptible[:,i] == sa.susceptible[:,i-1]
      sa.event_type[i] = "R"
    end
    sa.new_infection[:,i] = sa.infectious[:,i] & sa.susceptible[:,i-1]
    if size(event_db.events)[2]==4
      sa.new_recovery[:,i] = sa.recovered[:,i] & sa.infectious[:,i-1]
    end
  end
  return sa
end


"""
Determine recovery times for individuals
"""
function find_recovery_times(event_db::edb, narm=true)
  if narm == false
    return event_db.events[:,4]  - event_db.events[:,3]
  end
  if narm == true
    times = event_db.events[:,4]  - event_db.events[:,3]
    return times[isnan(times).==false]
  end
end

function create_loglikelihood(pop_db::DataFrame, event_db::edb)
  """
  Create an efficient log likelihood function for continuous and discrete
  SI and SIR models
  """
  sa = state_array(event_db)
  distance_mat = create_dist_mtx(pop_db)
  if event_db.cd == "continuous" && size(event_db.events)[2]==4
    """
    loglikleihood for continuous SIR model
    """
    function cSIR_ll(parameters)
      if parameters[1] <= 0 || parameters[2] <= 0 || parameters[3] <= 0
        return -Inf
      else
        dist_ab = dist_ab_mtx(distance_mat, parameters[1], parameters[2])
        ilm_ll = 0.0
        for i = 2:length(sa.event_times)
          if sum(sa.new_infection[:,i]) > 0
            exp_rates=sum(dist_ab[sa.new_infection[:,i], sa.infectious[:,i-1]], 2).^-1.0
            for j = 1:length(exp_rates)
              ilm_ll += loglikelihood(Exponential(exp_rates[j]), [sa.event_times[i]-sa.event_times[i-1]])
            end
          end
          if sum(sa.new_recovery[:,i]) > 0
            ilm_ll += loglikelihood(Exponential(parameters[3]), fill(sa.event_times[i]-sa.event_times[i-1], sum(sa.new_recovery[:,i])))
          end
          ilm_ll += sum(dist_ab[sa.susceptible[:,i], sa.infectious[:,i-1]])*(sa.event_times[i-1]-sa.event_times[i])
          ilm_ll += sum(sa.infectious[:,i]) * (1/parameters[3]) * (sa.event_times[i-1]-sa.event_times[i])
        end
        return ilm_ll
      end
    end
    return cSIR_ll
  end

  if event_db.cd == "continuous" && size(event_db.events)[2]==3
    """
    loglikleihood for continuous SI model
    """
    function cSI_ll(parameters)
      if parameters[1] <= 0 || parameters[2] <= 0
        return -Inf
      else
        dist_ab = dist_ab_mtx(distance_mat, parameters[1], parameters[2])
        ilm_ll = 0.0
        for i = 2:length(sa.event_times)
            exp_rates=sum(dist_ab[sa.new_infection[:,i], sa.infectious[:,i-1]], 2).^-1.0
            for j = 1:length(exp_rates)
              ilm_ll += loglikelihood(Exponential(exp_rates[j]), [sa.event_times[i]-sa.event_times[i-1]])
            end
          ilm_ll += sum(dist_ab[sa.susceptible[:,i], sa.infectious[:,i-1]])*(sa.event_times[i-1]-sa.event_times[i])
        end
        return ilm_ll
      end
    end
    return cSI_ll
  end

  if event_db.cd == "discrete" && size(event_db.events)[2]==4
    rt = find_recovery_times(event_db)
    """
    loglikleihood for discrete SIR model
    """
    function dSIR_ll(parameters)
      if parameters[1] <= 0 || parameters[2] <= 0 || parameters[3] <= 0
        return -Inf
      else
        dist_ab = dist_ab_mtx(distance_mat, parameters[1], parameters[2])
        ilm_ll = 0.0
        for i = 2:maximum(event_db.events[:,3])
          ilm_ll += -sum(dist_ab[sa.susceptible[:,i], sa.infectious[:,i-1]])
          ilm_ll += sum(log(1-exp(-sum(dist_ab[sa.new_infection[:,i], sa.infectious[:,i-1]], 2))))
        end
        ilm_ll += loglikelihood(Geometric(1/parameters[3]), rt-1)
        return ilm_ll
      end
    end
    return dSIR_ll
  end

  if event_db.cd == "discrete" && size(event_db.events)[2]==3
    """
    loglikleihood for discrete SI model
    """
    function dSI_ll(parameters)
      if parameters[1] <= 0 || parameters[2] <= 0
        return -Inf
      else
        dist_ab = dist_ab_mtx(distance_mat, parameters[1], parameters[2])
        ilm_ll = 0.0
        for i = 2:maximum(event_db.events[:,3])
          ilm_ll += -sum(dist_ab[sa.susceptible[:,i], sa.infectious[:,i-1]])
          ilm_ll += sum(log(1-exp(-sum(dist_ab[sa.new_infection[:,i], sa.infectious[:,i-1]], 2))))
        end
        return ilm_ll
      end
    end
    return dSI_ll
  end
end
