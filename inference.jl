type sdb
  susceptible::Array{Bool,2}
  infectious::Array{Bool,2}
  new_infection::Array{Bool,2}
  event_times::Array{Float64,1}
  event_type::Array{ASCIIString,1}
end

function state_array(event_db::edb)
  """
  Save the likelihood function from repetively determining state,
  by producing an array which contains this information for all 
  relevant time steps
  """
  sa=sdb(fill(false, (size(event_db.events)[1], length(unique(event_db.event_times[event_db.event_times.<Inf])))), fill(false, (size(event_db.events)[1], length(unique(event_db.event_times[event_db.event_times.<Inf])))), fill(false, (size(event_db.events)[1], length(unique(event_db.event_times[event_db.event_times.<Inf])))), unique(event_db.event_times[event_db.event_times.<Inf]), fill("I", length(unique(event_db.event_times[event_db.event_times.<Inf]))))
  for i = 1:length(sa.event_times)
    sa.susceptible[:,i] = find_state(event_db, sa.event_times[i], "S")
    sa.infectious[:,i] = find_state(event_db, sa.event_times[i], "I")
  end
  for i = 2:length(sa.event_times)
    if sa.susceptible[:,i] == sa.susceptible[:,i-1]
      sa.event_type[i] = "R"
    end
    sa.new_infection[:,i] = sa.infectious[:,i] .== sa.susceptible[:,i-1]
  end
  return sa
end

function find_recovery_times(event_db::edb, narm=true)
  """
  Determine recovery times for individuals
  """
  if narm == false
    return event_db.events[:,4]  - event_db.events[:,3]
  end
  if narm == true
    times = event_db.events[:,4]  - event_db.events[:,3]
    return times[isnan(times).==false]
  end
end

function create_loglikelihood(pop_db, event_db::edb)
  """
  Create an efficient log likelihood function for continuous and discrete 
  SI and SIR models
  """
  sa = state_array(event_db)
  distance_mat = create_dist_mtx(pop_db)
  if size(event_db.events)[2]==4
    rt = find_recovery_times(event_db)
  end
  if event_db.cd == "continuous" && size(event_db.events)[2]==4
    function cSIR_ll(parameters)
      """
      loglikleihood for continuous SIR model
      """
      if parameters[1] <= 0 || parameters[2] <= 0 || parameters[3] <= 0
        return -Inf
      else
        dist_ab = dist_ab_mtx(distance_mat, parameters[1], parameters[2])
        ilm_ll = loglikelihood(Exponential(parameters[3]), rt)
        for i = 2:length(sa.event_times)
          if sa.event_type == "I"
            exp_rates=sum(dist_ab[sa.new_infection[:,i], sa.infectious[:,i-1]], 2).^-1.0
            for j = 1:length(exp_rates)
              ilm_ll += loglikelihood(Exponential(exp_rates[j]), [sa.event_times[i]-sa.event_times[i-1]])
            end
            ilm_ll += sum(dist_ab[sa.susceptible[:,i], sa.infectious[:,i-1]])*(sa.event_times[i-1]-sa.event_times[i])
          end
        end
        return ilm_ll
      end
    end
    return cSIR_ll
  end

  if event_db.cd == "continuous" && size(event_db.events)[2]==3
    function cSI_ll(parameters)
      """
      loglikleihood for continuous SI model
      """
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
    function dSIR_ll(parameters)
      """
      loglikleihood for discrete SIR model
      """
      if parameters[1] <= 0 || parameters[2] <= 0 || parameters[3] <= 0
        return -Inf
      else
        dist_ab = dist_ab_mtx(distance_mat, parameters[1], parameters[2])
        ilm_ll = loglikelihood(Geometric(1/parameters[3]), rt)
        for i = 2:maximum(event_db.events[:,2])
        end
        return ilm_ll
      end
    end
    return dSIR_ll
  end

  if event_db.cd == "discrete" && size(event_db.events)[2]==3
    function dSI_ll(parameters)
      """
      loglikleihood for discrete SIR model
      """
      if parameters[1] <= 0 || parameters[2] <= 0
        return -Inf
      else
        dist_ab = dist_ab_mtx(distance_mat, parameters[1], parameters[2])
        ilm_ll = 0.0
        for i = 2:maximum(event_db.events[:,2])
        end
        return ilm_ll
      end
    end
    return dSI_ll
  end
end

#=
function create_sir_loglikelihood(distance_mat, susceptible_array, infectious_array, recovery_times)
"""
Create the SIR loglikelihood function, which accepts only 3 parameters,
alpha, beta, and gamma_inverse
"""
  function sir_loglikelihood(parameters)
  """
  Compute the log likelihood of discrete time SIR models (alpha, beta, gamma)
  """
    distance_mat_alphabeta = distance_mat_alphabeta_fun(distance_mat, parameters[1], parameters[2])
    if parameters[1] <= 0 || parameters[2] <= 0 || parameters[3] >= 1 || parameters[3] <= 0
      -Inf
    else
      daily_loglikes = zeros(sum(sum(susceptible_array,1) .> 0)-1)
      for i = 1:length(daily_loglikes)
        infect_probs=nans(size(susceptible_array)[1])
        infect_probs[susceptible_array[:,i]]=infect_prob_fun(distance_mat_alphabeta, infectious_array[:,i], susceptible_array[:,i])
        daily_loglikes[i]=sum(log([(infect_probs[susceptible_array[:,i] & infectious_array[:,i+1]]), 1.-(infect_probs[susceptible_array[:,i] & susceptible_array[:,i+1]])]))
      end
      loglikelihood(Geometric(parameters[3]), recovery_times) + logpdf(Uniform(0,1), parameters[3]) + logpdf(Gamma(1,1), parameters[1]) + logpdf(Gamma(5,1), parameters[2]) + sum(daily_loglikes)
    end
  end
sir_loglikelihood
end
=#