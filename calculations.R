library(deSolve)
library(cubature)

expectation <- function(fn) {
  expectation <- integrate(
    \(x) x * fn(x),
    0,
    50,
    abs.tol = 1e-4,
    rel.tol = 1e-4
  )
  
  expectation$value
}

multi_forward_gt <- function(m, strain) {
  # effective S
  fn.S <- \(t) m$S.AB(t) + x * m[[paste("S.", strain, sep="")]](t)
  fn.beta <- \(t) m$params[[strain]](0)$beta
  fn.g0 <- \(tau) m[[paste("g.0.", strain, sep="")]](tau)
  
  function(t, tau) {
    num <- fn.g0(tau) * fn.beta(t + tau) * fn.S(t + tau)
    denom <- integrate(
      \(sigma) fn.g0(sigma) * fn.beta(t + sigma) * fn.S(t + sigma),
      0,
      Inf,
      abs.tol = 1e-8
    )$value
    
    if (denom <= 0) {
      return(0)
    }
    
    num/denom
  }
}

multi_backward_gt <- function(m, strain) {
  # effective S
  fn.i <- function(t) {
    result <- m[[paste("incidence.", strain, "1", sep="")]](t) + m[[paste("incidence.", strain, "2", sep="")]](t)
    # print("START")
    # print(t)
    # print(result)
  }
  fn.beta <- \(t) m$params[[strain]](0)$beta
  fn.g0 <- \(tau) m[[paste("g.0.", strain, sep="")]](tau)
  
  function(t, tau) {
    num <- fn.g0(tau) * fn.i(t-tau)
    
    denom <- integrate(
      \(sigma) fn.g0(sigma) * fn.i(t-sigma),
      0,
      100,
      abs.tol = 1e-4,
      rel.tol = 1e-4
    )$value
    
    if (denom <= 0) {
      return(0)
    }
    
    num/denom
    
  }
}

# return a function which calculates the intrinsic generation time for a single
# model strain
multi_intrinsic_gt <- function(m, strain) {
  function(tau) {
    # this is just a hypoexponential distribution
    sigma.A <- m$params[[strain]](0)[["sigma"]]
    gamma.A <- m$params[[strain]](0)[["gamma"]]
    
    (sigma.A * gamma.A) * (exp(-gamma.A * tau) - exp(-sigma.A * tau)) / (sigma.A - gamma.A)
  }
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