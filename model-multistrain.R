library(deSolve)

source("solver.R")
source("simulation.R")

individual_diff_fn <- function(t, s, p) {
  with(as.list(c(s, p)), {
    p <- p(t)
    
    dL <-  - p[["sigma"]] * s[["L"]]
    dF <-  p[["sigma"]] * s[["L"]] - p[["gamma"]] * s[["F"]]
    list(c(dL, dF))
  })
}

seir_2_diff_fn <- function(t, state, parameters) {
  with(as.list(c(t, state, parameters)), {
    A <- A(t)
    B <- B(t)
    
    states = names(state)
    
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

simulate_two_seir <- function(initial_value, params, start, end) {
  # use a deterministic model
  data <- data.frame(ode(
    y = initial_value,
    func = seir_2_diff_fn,
    times = seq(start, end, by=1),
    parms = parameters,
    events = list(data = params$events)
  ))
  
  if (is.null(params$events) || nrow(params$events) == 0) {
    # empty - add names just to prevent issues
    params$events <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(params$events) <- c("var", "time", "value", "method")
  }
  
  result <- generate_summary(data, params)
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

generate_summary <- function(data, params) {
  result <- list()
  
  result$data <- data
  result$params <- params
  result$events <- params$events
  
  A_data <- data.frame(
    time = data$time,
    S = data$S.AB + data$S.A,
    E = data$E.A1 + data$E.A2,
    I = data$I.A1 + data$I.A2,
    R = data$R.A + data$S.B + data$E.B2 + data$I.B2 + data$R.AB,
    strain = "A"
  )
  
  B_data <- data.frame(
    time = data$time,
    S = data$S.AB + data$S.B,
    E = data$E.B1 + data$E.B2,
    I = data$I.B1 + data$I.B2,
    R = data$R.B + data$S.A + data$E.A2 + data$I.A2 + data$R.AB,
    strain = "B"
  )
  
  overview <- rbind(A_data, B_data)
  
  strains <- list()
  
  for (s in c("A", "B")) {
    strain <- list()
    strain$overview <- calculate_incidence(overview[overview$strain == s,])
    strain$params <- \(t) c(params$A(t), events=params$events, N=params$N)
    strain <- generate_approximations(strain)
    
    strains[[s]] = strain
  }
  
  result$overview <- overview
  result$strains <- strains
  
  return(result)
}

generate_approximations <- function(data) {
  data$individual <- data.frame(ode(
    y = c("L" = 1, "F" = 0),
    times = seq(0, 100, by = 0.1),
    parms = data$params,
    func = individual_diff_fn
  ))
  
  # set up approximation functions for S, E, I, R, incidence
  data$S <- approxfun(data$overview$time, data$overview$S, rule=2)
  data$E <- approxfun(data$overview$time, data$overview$E, rule=2)
  data$I <- approxfun(data$overview$time, data$overview$I, rule=2)
  data$R <- approxfun(data$overview$time, data$overview$R, rule=2)
  data$i <- approxfun(data$overview$time, data$overview$incidence, rule=2)
  
  data$F <- approxfun(data$individual$time, data$individual$F, yright=0, rule=2)
  data$L <- approxfun(data$individual$time, data$individual$L, yright=0, rule=2)
  
  return(data)
}