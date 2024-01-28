library(deSolve)

expectation <- function(fn) {
  expectation <- integrate(
    \(x) x * fn(x),
    0,
    Inf
  )
  
  expectation$value
}

# return a function which calculates the intrinsic generation time for a single
# model strain
single_intrinsic_gt <- function(m) {
  function(tau) {
    # this is just a hypoexponential distribution
    sigma <- m$params[["sigma"]]
    gamma <- m$params[["gamma"]]
    
    (sigma * gamma) * (exp(-gamma * tau) - exp(-sigma * tau)) / (sigma - gamma)
  }
}

# return a function which calculates the forward generation time for a single
# model strain
single_forward_gt <- function(m) {
  function(t, tau) {
    num <- m$g.0(tau) * m$beta(t + tau) * m$S(t + tau)
    denom <- integrate(
      \(sigma) m$g.0(sigma) * m$beta(t + sigma) * m$S(t + sigma),
      0,
      Inf
    )$value
    
    if (denom <= 0) {
      return(0)
    }
    
    num/denom
  }
}

single_backward_gt <- function(m) {
  function(t, tau) {
    num <- m$g.0(tau) * m$i(t - tau)
    denom <- integrate(
      \(sigma) m$g.0(sigma) * m$i(t - sigma),
      0,
      t
    )$value
    
    if (denom <= 0) {
      return(0)
    }
    
    num/denom
  }
}