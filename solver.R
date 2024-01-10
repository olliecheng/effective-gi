solve_stoch_de <- function(func, rng_func, initial_value, params, start, end, simulation=NULL) {
  state <- initial_value
  vars <- length(initial_value)
  
  result <- data.frame(time = start, as.list(initial_value))
  
  simulate <- !is.null(simulation)
  
  times = seq(start + 1, end, by = 1)
  for (t in times) {
    iteration <- func(state, params)
    
    for (i in seq_along(iteration)) {
      # calculate the difference randomly, and bound by the size of the
      # current state
      difference <- min(c(rng_func(iteration[[i]]), state[[i]]))
      
      # subtract from the current state
      state[[i]] <- state[[i]] - difference;
      
      # add to the next state
      state[[i+1]] <- state[[i+1]] + difference;
      
      # update the population
      if (simulate) {
        if (difference != 0) {
          simulation$data <- simulation$transitions[[i]](
            simulation$data,
            difference,
            t
          )
        }
        
        results <- sapply(c(1,2,3,4), \(x) sum(simulation$data$.state == x))
        
        if(results[[1]] != state[[1]] || results[[2]] != state[[2]] ||
           results[[3]] != state[[3]] || results[[4]] != state[[4]]) {
          print(i)
          print(difference)
          
          print("IMBALANCE")
          
          print(results)
          print(state)
          flkasdfldsakjfads
        }
      }
    }
    
    # add result to list
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
    rv$pop <- simulation$data
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