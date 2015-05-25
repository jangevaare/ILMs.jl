module ilmtools

using DataFrames, Distances, Distributions, Gadfly, Interact, Mamba

export
  create_event_db,
  create_pop_db_univariates,
  create_pop_db_bivariate,
  infect_recover_loop,
  state_timeseries,
  state_animation,
  create_loglikelihood

include("population.jl")
include("events.jl")
include("visuals.jl")
include("inference.jl")

end # module
