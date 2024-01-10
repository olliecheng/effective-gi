initialise_simulation <- function(initial_value, parameters, transition_fn) {
  # create a matrix of the population
  initial <- matrix(0, parameters[["N"]], length(initial_value))
  population <- data.frame(initial)
  colnames(population) <- names(initial_value) # set column names
  
  population$.infected_by = NA
  population$.infected = NA
  
  # give each member a default starting state
  population$.state <- do.call(
    c,
    sapply(
      seq_along(initial_value),
      function(v) replicate(initial_value[[v]], v),
      simplify = FALSE
    )
  )
  
  return <- c()
  return$data <- population
  return$transitions <- transition_fn
  
  return
}

initialise_seir_simulation <- function(initial_value, parameters) {
  initialise_simulation(initial_value, parameters, seir_sim_transitions)
}

seir_sim_transitions <- c(
  # infected for the first time, generate forwards and backwards interval
  S_to_E = function(population, difference, t) {
    subset = which(population$.state == 1)
    subset_sample = sample(subset, difference, replace=FALSE)
    
    stopifnot(nrow(subset_sample) == difference)
    
    population[subset_sample, ][,2] <- t
    population[subset_sample, ]$.state <- 2
    
    population
  },
  E_to_I = function(population, difference, t) {
    subset = which(population$.state == 2)
    subset_sample = sample(subset, difference, replace=FALSE)
    
    stopifnot(nrow(subset_sample) == difference)
    
    population[subset_sample, ][,3] <- t
    population[subset_sample, ]$.state <- 3
    
    population
  },
  I_to_R = function(population, difference, t) {
    subset = which(population$.state == 3)
    subset_sample = sample(subset, difference, replace=FALSE)
    
    stopifnot(nrow(subset_sample) == difference)
    
    population[subset_sample, ][,4] <- t
    population[subset_sample, ]$.state <- 4
    
    population
  }
)

get_counts <- function(population) {
  sum(u, na.rm=TRUE)
}