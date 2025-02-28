nrep = 2 ## num reps
out = matrix(NA, ncol = 6, nrow = nrep) #store optimization

for(r in 1:nrep){
  
  print(r) #displays iteration
  
  # set parameters
  P = 365
  d = 10
  n = P * d
  
  mu_C = -0.75
  
  beta_0 = 1 #amplitude of seasonal effect
  beta_1 = 0
  tau = 150
  
  
  beta_2 = 1
  beta_3 = .5
  tau_C = 120
  Seasonal_C = beta_2 +beta_3 * cos((2*pi*(1:n - tau_C))/P)
  plot(Seasonal_C)
  
  sig_C = 1
  sig_X = 1
  sig_TT = 25
  
  e_TT = rnorm(n, sd = sig_TT)
  TT = 32 - 25 * cos((2*pi*(1:n - 120))/P) + e_TT |>
    ts(start = c(1990, 1), frequency = P)
  
  # Draw changes
  eps_C = rnorm(n, mean = 0, sd = Seasonal_C)
  C = mu_C + beta_0 * cos((2*pi*(1:n - tau))/P) + beta_1 * TT + eps_C |>
    ts(start = c(1990, 1), frequency = P)
  
  plot(C)
  abline(a = mu_C, b = 0)
  
  
  # Create Xt as function of changes
  X = rep(0, n)
  for(t in 2:n) X[t] = max(X[t-1] + C[t], 0 )
  X = ts(X, start = c(1990, 1), frequency = P)
  
  if(F){
    op = par(mfrow = c(2, 1))
    plot(X, ylim = range(X, C)); abline(h = 0, lty = 'dotted')
    plot(C, ylim = range(X, C)); abline(h = 0, lty = 'dotted')
    par(op)
  }
  
  
  dC = function(x, theta, t){
    my_mu = m(theta, t = t)
    my_sig = theta[4] +theta[5] * cos((2*pi*(t - theta[6]))/P)
    dnorm(x = x, mean = my_mu, sd = my_sig)
  }
  
  pC = function(q, theta, t){
    my_mu = m(theta, t = t)
    my_sig = theta[4] +theta[5] * cos((2*pi*(t - theta[6]))/P)
    pnorm(q = q, mean = my_mu, sd = my_sig)
  }
  
  lik = function(theta, X){
    out = 0
    for(t in 2:n){
      if(X[t] > 0){
        new_part = dC(x = X[t] - X[t-1], theta, t)
      } else {
        new_part = pC(q = -X[t-1], theta, t)
      }
      out = out + log(new_part)
    }
    return(-1*out)
  }
  
  
  # # Plot lik
  # if(FALSE){
  #   xax = seq(-2, 2, .1)
  #   f = rep(NA, length(xax))
  #   for(i in seq_along(xax)){
  #     print(xax[i])
  #     f[i] = lik(theta = c(mu_C, xax[i], tau, 1), X = X)
  #   }
  #   plot(f ~ xax)
  # }
  
  optim_out = optim(c(-.75, 1, 150, 1,0.5,120),
                    fn = lik,
                    X = X,
                    lower = c(-1.5, .9, 100, .8,.3,119),
                    upper = c(-.5, 1.1, 200, 1.1,.7,121),
                    method = "L-BFGS-B")
  out[r, ] = optim_out$par
}

theta_true = c(-.75, 1, 150, 1,.5,120)
names = c("mu", "beta0", "tau", "beta2","beta3","tauc")
op = par(mfrow = c(3, 2))
for(i in 1:6){
  boxplot(out[, i], main = names[i])
  abline(h = theta_true[i])
}
par(op)

####2

library(tsbox)

out = replicate(n = 30, expr = {
  
  # set parameters
  P = 365
  d = 100
  n = P * d
  
  mu_C = -0.75
  
  beta_0 = 0
  beta_1 = 0
  tau = 150
  
  sig_C = 1
  sig_X = 1
  sig_TT = 25
  
  e_TT = rnorm(n, sd = sig_TT)
  TT = 32 - 25 * cos((2*pi*(1:n - 120))/P) + e_TT |>
    ts(start = c(1990, 1), frequency = P)
  
  # Draw changes
  eps_C = rnorm(n, mean = 0, sd = sig_C)
  C = mu_C + beta_0 * cos((2*pi*(1:n - tau))/P) + beta_1 * TT + eps_C |>
    ts(start = c(1990, 1), frequency = P)
  
  # Create Xt as function of changes
  X = rep(0, n)
  for(t in 2:n) X[t] = max(X[t-1] + C[t], 0 )
  X = ts(X, start = c(1990, 1), frequency = P)
  
  C_est = c(0, diff(X))
  
  if(FALSE){
    op = par(mfrow = c(3, 1))
    plot(X, ylim = range(X, C, C_est)); abline(h = 0, lty = 'dotted')
    plot(C, ylim = range(X, C, C_est)); abline(h = 0, lty = 'dotted')
    plot(C_est, ylim = range(X, C, C_est)); abline(h = 0, lty = 'dotted')
    par(op)
  }
  
  # Estimate the mean and sd of the C process using the estimated C = diff(X)
  C_est_mean = matrix(C_est, nrow = P) |> rowMeans()
  C_est_sd = matrix(C_est, nrow = P) |> apply(1, sd)
  
  
  if(F){
    op = par(mfrow = c(2, 1))
    plot(C_est_sd, ylim = range(0, sig_C, C_est_sd),
         main = "Estimated SD from hat{C} = diff(X)")
    abline(h = sig_C, col = 'red', lty = 'dotted')
    abline(h = 0, lty = 'dotted')
    
    C_sd = matrix(C, nrow = P) |> apply(1, sd)
    plot(C_sd, ylim = range(0, sig_C, C_sd),
         main = "Estimated SD from true C (unobservable)")
    abline(h = sig_C, col = 'red', lty = 'dotted')
    abline(h = 0, lty = 'dotted')
    par(op)
  }
  
  # sum or squares as function of mu
  # Equation (4) from paper
  SS = function(theta) {
    S = 0
    for (t in 2:n) {
      # sig = max(C_est_sd[t2nu(t)], .05)
      sig = 1
      EE = E(X[t - 1], theta, t, m, sig)
      S = S + (X[t] - EE) ^ 2 / sig ^ 2
    }
    return(S)
  }
  
  # optimization
  # optim_fit = optim(-.25, SS, method = "Brent", lower = -2, upper = 0)
  # cbind(estimated = optim_fit$par, true = mu_C)
  
  # get parameter estimates
  mu_hat = optim_fit$par
  
  # Xhat = E[X_t | X_{t-1}]
  Xhat = rep(NA, n)
  Xhat[1] = 0
  for(t in 2:n){
    Xhat[t] = E(X = X[t - 1], theta = mu_C, t = t, m = m, sig_C = sig_C)
  }
  
  # get residuals
  r = X - Xhat
  
  if(F) {
    op = par(mfrow = c(3, 1))
    plot.ts(X[1:365], main = "True X")
    plot.ts(Xhat[1:365], main = "hat{X}")
    plot.ts(r[1:365], main = "residual")
    par(op)
  }
  
  # Try to write sum-of-squares to estimates sig_C as
  # each term is rt - E[Xt - mu | X_{t-1}]
  SS_sig = function(sig_C){
    
    m = mu_C
    
    S = 0
    for(t in 2:n){
      
      k =  (X[t-1] + m) / sig_X
      
      A = E2(X[t-1], mu_C, sig_C)
      
      B = E(X[t - 1], mu_C, t, m, sig_C)
      
      D = A - B^2
      
      S = S + (r[t-1]^2 - D)^2
    }
    return(S)
  }
  
  # Plot SS_sig
  if (F) {
    xax = seq(0, 2, .2)
    f = rep(NA, length(xax))
    for (i in seq_along(xax)) {
      print(xax[i])
      f[i] = SS_sig(xax[i])
    }
    plot(f ~ xax)
  }
  
  # optimization
  optim_fit2 = optim(1, SS_sig, method = "Brent", lower = 0, upper = 5)
  # cbind(estimated = optim_fit2$par, true = sig_C)
  optim_fit2$par
  
 ###### 3
  
   
  t2nu = function(t, P = 365){
    ((t-1) %% P) + 1
  }
  
  
  t2d = function(t, P = 365){
    floor((t-1)/P)
  }
  
  dnu2t = function(d, nu, P = 365){
    if(nu > P) stop("nu exceeds P")
    d*P + nu
  }
  
  
  
  ### additional 
  m = function(theta, t, P = 365) {
    
    nu = t2nu(t)
    
    # extract parameters from theta for readability
    mu = theta[1]
    beta_0 = theta[2]
    tau = theta[3]
    
    # return mean function
    mu + beta_0 * cos(2 * pi * (nu - tau) / P)
  }
  
