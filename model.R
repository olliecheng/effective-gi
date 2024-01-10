source("solver.R")

seir_diff_fn <- function(state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <-  -( (beta / N) * S * I )
    dE <- (beta / N) * S * I - sigma * E
    dI <- -( gamma * I ) + sigma * E
    dR <- gamma * I
    c(dS, dE, dI, dR)
  })
}

seir_stoch_diff_fn <- function(state, parameters) {
  with(as.list(c(state, parameters)), {
    S_to_E <- (beta / N) * S * I
    E_to_I <- sigma * E
    I_to_R <- gamma * I
    c(S_to_E, E_to_I, I_to_R)
  })
}

simulate_seir <- function(initial_value, params, start, end, stochastic = c(TRUE, FALSE), simulate=c(TRUE, FALSE)) {
  # by default, simulate a stochastic model
  if (missing(stochastic) || stochastic) {
    result <- solve_stoch_de(
        seir_stoch_diff_fn,
        \(x) rpois(1, x), # randomly generate from a Poisson distribution
        initial_value = state,
        params = parameters,
        start = start,
        end = end,
        simulate = missing(simulate) || simulate
      )
  } else {
    result <- solve_de(
      seir_diff_fn,
      initial_value = state,
      params = parameters,
      start = start,
      end = end,
      dt = 0.1
    )
    
    # https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
    is.wholenumber <- function(x, tol = 1e-5) {
      min(abs(c(x%%1, x%%1-1))) < tol
    }
    
    whole <- sapply(result$overview$time, is.wholenumber)
    result$overview <- result$overview[whole, ]
  }
  
  result$overview <- calculate_incidence(result$overview)
  
  return(result)
}


calculate_incidence <- function(df) {
  df$incidence <- c(0, sapply(diff(df$S), \(x) -x))
  df
}