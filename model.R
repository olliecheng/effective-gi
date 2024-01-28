library(deSolve)

source("solver.R")
source("simulation.R")

individual_diff_fn <- function(t, s, p) {
  with(as.list(c(s, p)), {
    dL <-  - p[["sigma"]] * s[["L"]]
    dF <-  p[["sigma"]] * s[["L"]] - p[["gamma"]] * s[["F"]]
    list(c(dL, dF))
  })
}

seir_diff_fn <- function(t, state, parameters) {
  with(as.list(c(t, state, parameters)), {
    beta <- beta_fn(t)
    
    dS <-  -( (beta / N) * S * I )
    dE <- (beta / N) * S * I - sigma * E
    dI <- -( gamma * I ) + sigma * E
    dR <- gamma * I
    list(c(dS, dE, dI, dR))
  })
}

simulate_seir <- function(initial_value, params, start, end) {
  # turn beta_events into a format recognised by deSolve
  params$beta_fn <- approxfun(params$beta_events$time, params$beta_events$value, method="constant", rule=2)
  
  result <- list()
  result$overview <- data.frame(ode(
    y = initial_value,
    func = seir_diff_fn,
    times = seq(start, end, by=0.1),
    parms = params,
  ))
  
  result$individual <- data.frame(ode(
    y = c("L" = 1, "F" = 0),
    times = seq(0, 50, by = int),
    parms = parameters,
    func = individual_diff_fn
  ))
  
  
  result$overview <- calculate_incidence(result$overview)
  result$params <- parameters
  
  result$beta <- params$beta_fn
  result$events <- tail(params$beta_events, -1) # everything except the first
  
  # set up approximation functions for S, E, I, R, incidence
  result$S <- approxfun(result$overview$time, result$overview$S, yright=0, rule=2)
  result$E <- approxfun(result$overview$time, result$overview$E, rule=2)
  result$I <- approxfun(result$overview$time, result$overview$I, rule=2)
  result$R <- approxfun(result$overview$time, result$overview$R, rule=2)
  result$i <- approxfun(result$overview$time, result$overview$incidence, rule=2)
  
  result$F <- approxfun(result$individual$time, result$individual$F, yright=0, rule=2)
  result$L <- approxfun(result$individual$time, result$individual$L, yright=0, rule=2)
  
  return(result)
}

calculate_incidence <- function(df) {
  df$incidence <- c(
    0, # starting value
    sapply(
      diff(df$S) / diff(df$time),
      \(x) -x
    )
  )
  
  df
}