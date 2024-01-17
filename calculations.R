library(deSolve)

basic_mean_gt <- function(m) {
  num <- integrate(
    \(tau) tau * m$F(tau),
    0, Inf
  )
  
  denom <- integrate(
    \(tau) m$F(tau),
    0, Inf
  )
  
  if (denom$value == 0) {
    return(0)
  }
  
  num$value / denom$value
}

forward_mean_gt <- function(m, exposure) {
  num <- integrate(
      \(tau) tau * m$F(tau) * m$S(exposure + tau),
      0,
      Inf
    )
  
  denom <- integrate(
    \(tau) m$F(tau) * m$S(exposure + tau),
    0,
    Inf
  )
  
  if (denom$value == 0) {
    return(0)
  }
  
  num$value / denom$value
}

backward_mean_gt <- function(m, exposure) {
  num <- integrate(
    \(tau) tau * m$F(tau) * m$i(exposure - tau),
    0,
    Inf
  )
  
  denom <- integrate(
    \(tau) m$F(tau) * m$i(exposure - tau),
    0,
    Inf
  )
  
  if (denom$value == 0) {
    return(0)
  }
  
  num$value / denom$value
}