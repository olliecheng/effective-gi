library(deSolve)

source("solver.R")
source("simulation.R")

eventfun <- function(t, y, parms){
  with(as.list(y), {
    y[y < 1e-6] <- 0
    return(y)
  })
}

seir_2_diff_fn <- function(t, state, parameters) {
  with(as.list(c(t, state, parameters)), {
    # format:
    # <strain> . <stage> . <start_state> . <end_state>
    
    states = c("S.AB", "S.A", "S.B", "E.A1", "E.B1", "E.A2", "E.B2",
               "I.A1", "I.B1", "I.A2", "I.B2", "R.A", "R.B", "R.AB")
    
    # matrix of all the transitions, from x to y is accessed by t[x, y]
    t <- matrix(nrow = 14, ncol = 14, dimnames=list(from=states, to=states))
    
    # susceptible -> exposed
    t["S.AB", "E.A1"] <-     (A$beta / N) * S.AB * (I.A1 + I.A2)
    t["S.AB", "E.B1"] <-     (B$beta / N) * S.AB * (I.B1 + I.B2)
    t["S.A",  "E.A2"] <- x * (A$beta / N) * S.A  * (I.A1 + I.A2)
    t["S.B",  "E.B2"] <- x * (B$beta / N) * S.B  * (I.B1 + I.B2)
    
    # exposed -> infectious
    t["E.A1", "I.A1"] <- A$sigma * E.A1
    t["E.A2", "I.A2"] <- A$sigma * E.A2
    t["E.B1", "I.B1"] <- B$sigma * E.B1
    t["E.B2", "I.B2"] <- B$sigma * E.B2
    
    # infectious -> recovered
    t["I.A1", "R.A" ] <- A$gamma * I.A1
    t["I.A2", "R.AB"] <- A$gamma * I.A2
    t["I.B1", "R.B" ] <- B$gamma * I.B1
    t["I.B2", "R.AB"] <- B$gamma * I.B2
    
    # recovered from A -> susceptible to B
    #  i.e. temporary complete immunity period
    t["R.A",  "S.B" ] <- A$omega * R.A
    t["R.B",  "S.A" ] <- B$omega * R.B
    
    # sum up all the transitions from the matrix t
    derivatives <- list(
      sapply(
        states, function(x) {
          # get everything LEAVING this state
          leaving  <- sum(t[x, ], na.rm=TRUE)
          entering <- sum(t[, x], na.rm=TRUE)
          
          # total rate is the difference
          return(entering - leaving)
        }
      )
    )
    
    return(derivatives)
  })
}

simulate_two_seir <- function(initial_value, params, start, end, stochastic=TRUE, simulate=TRUE) {
  # by default, simulate a stochastic model
  if (stochastic) {
  }
  
  # use a deterministic model
  else {
    result <- list()
    data <- data.frame(ode(
      y = initial_value,
      func = seir_2_diff_fn,
      times = seq(start, end, by=1),
      parms = parameters,
      #method = "euler",
      #atol = 1e-2,
      #events = list(func = eventfun, time = start:end)
    ))
    result$data <- data
    
    result$overview <- data.frame(
      time = data$time,
      S.A = data$S.AB + data$S.A,
      S.B = data$S.AB + data$S.B,
      E.A = data$E.A1 + data$E.A2,
      E.B = data$E.B1 + data$E.B2,
      I.A = data$I.A1 + data$I.A2,
      I.B = data$I.B1 + data$I.B2,
      R = data$R.AB
    )
    
    result$strains <- list()
    result$strains$A <- list()
    result$strains$A$overview <- data.frame(
      time = data$time,
      S = data$S.AB + data$S.A,
      E = data$E.A1 + data$E.A2,
      I = data$I.A1 + data$I.A2,
      R = data$R.A + data$S.B + data$E.B2 + data$I.B2 + data$R.AB
    )
    
    result$strains$B <- list()
    result$strains$B$overview <- data.frame(
      time = data$time,
      S = data$S.AB + data$S.B,
      E = data$E.B1 + data$E.B2,
      I = data$I.B1 + data$I.B2,
      R = data$R.B + data$S.A + data$E.A2 + data$I.A2 + data$R.AB
    )
    
    return(result)
    
    result$individual <- data.frame(ode(
      y = c("L" = 1, "F" = 0),
      times = seq(0, 50, by = int),
      parms = parameters,
      func = individual_diff_fn
    ))
  }
  
  result$overview <- calculate_incidence(result$overview)
  
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