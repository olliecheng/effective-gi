solve_stoch_de <- function(func, rng_func, initial_value, params, start, end, simulate=c(TRUE, FALSE)) {
  state <- initial_value
  vars <- length(initial_value)
  
  result <- data.frame(time = start, as.list(initial_value))
  
  # create a matrix of the population
  initial <- matrix(0, parameters[["N"]], vars)
  population <- data.frame(initial)
  colnames(population) <- names(state) # set column names
  
  # give each member a default starting state
  population$.state <- do.call(
    c,
    sapply(
      seq_along(state),
      function(v) replicate(state[[v]], v),
      simplify = FALSE
    )
  )
  
  times = seq(start + 1, end, by = 1)
  for (t in times) {
    iteration <- func(state, params)
    
    for (i in seq_along(iteration)) {
      # calculate the difference randomly, and bound by the size of the
      # current state
      difference <- min(c(rng_func(iteration[[i]]), state[[i]]))
      
      # subtract from the current index
      # ensure that the state is non-negative; 
      # it should never be, but this is just to check
      state[[i]] <- state[[i]] - difference;
      
      # add to the next state
      state[[i+1]] <- state[[i+1]] + difference;
      
      # update the population
      if (simulate) {
        if (difference > 0) {
          subset = which(population$.state == i)
          subset_sample = sample(subset, difference)
          
          stopifnot(nrow(subset_sample) == difference)
          
          population[subset_sample, ][,(i+1)] <- t
          population[subset_sample, ]$.state <- i+1
        }
      }
    }
    
    result <- rbind(
      result,
      as.list(
        c(state, c(time = t))
      )
    )
  }
  
  rv <- list()
  rv$overview <- result
  if (simulate) {
    rv$pop <- population
  } else {
    rv$pop <- NULL
  }
  
  return(rv)
}

# Function to solve a differential equation by the Euler Method
# The input function must return dX in the same index position as X in
# the initial_value vector. All params are passed into the input function.
solve_de <- function(func, initial_value, params, start, end, dt) {
  state <- initial_value
  result <- data.frame(time = start, as.list(initial_value))
  
  times = tail(seq(start, end, by = dt), n = -1)
  for (t in times) {
    iteration <- func(state, params) |> 
      sapply(\(x) x * dt) # multiply all elements by dt
    
    # apply Euler's Method i.e. X <- X + dX/dt
    state <- rowSums(cbind(state, iteration))
    
    result <- rbind(
      result,
      as.list(
        c(state, c(time = t))
      )
    )
  }
  
  rv <- list()
  rv$overview <- result
  rv$pop <- NULL # no population simulation has occurred here
  
  return(rv)
}