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
  return$data <- c()
  return$data$pop <- population
  return$data$events <- data.frame("infector"=integer(), "infectee"=integer(), "time"=integer())
  return$transitions <- transition_fn
  
  return
}

initialise_seir_simulation <- function(initial_value, parameters) {
  initialise_simulation(initial_value, parameters, seir_sim_transitions)
}

seir_sim_transitions <- c(
  # infected for the first time, generate forwards and backwards interval
  S_to_E = function(data, difference, t) {
    subset = which(data$pop$.state == 1)
    subset_sample = sample_vec(subset, difference, replace=FALSE)
    
    stopifnot(nrow(subset_sample) == difference)
    
    data$pop[subset_sample, ][,2] <- t
    data$pop[subset_sample, ]$.state <- 2
    
    # generate event information
    infectors <- sample_vec(
      which(data$pop$.state == 3),
      difference,
      replace=TRUE
    )
    events <- data.frame(infector = infectors, infectee = subset_sample, time = t)
    data$events <- rbind(data$events, events)
    
    data
  },
  
  E_to_I = function(data, difference, t) {
    subset = which(data$pop$.state == 2)
    subset_sample = sample_vec(subset, difference, replace=FALSE)
    
    stopifnot(nrow(subset_sample) == difference)
    
    data$pop[subset_sample, ][,3] <- t
    data$pop[subset_sample, ]$.state <- 3
    
    data
  },
  
  I_to_R = function(data, difference, t) {
    subset = which(data$pop$.state == 3)
    subset_sample = sample_vec(subset, difference, replace=FALSE)
    
    stopifnot(nrow(subset_sample) == difference)
    
    data$pop[subset_sample, ][,4] <- t
    data$pop[subset_sample, ]$.state <- 4
    
    data
  }
)

get_counts <- function(population) {
  sum(u, na.rm=TRUE)
}

sample_vec <- function(v, n, replace=FALSE) {
  # if there's only one thing to return, return it
  if (length(v) == 1) {
    # check that they're not looking for >1 replications,
    # unless replace is TRUE
    
    stopifnot(n == 1 || replace)
    return(replicate(n, v))
  } else {
    # return a random sample
    return(sample(v, n, replace))
  }
}