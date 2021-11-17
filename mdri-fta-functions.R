
# MDRI-FTA

# Primary function -----------------------------------------------------------

f.est.adjmdri <- function(bigt 
                          , mdri.unadj 
                          , mdri.shape 
                          , p.vs.fixed 
                          , vs.fixed
                          , p.vs.var  
                          , vs.x1, vs.p1 
                          , vs.x2, vs.p2 
                          , p.i 
                          , delaytotest.i 
                          , delta.mean.i 
                          , p.r 
                          , delta.fixed.r 
                          , p.f 
                          , delta.min.f 
                          , delta.max.f 
                          , p.n 
                          , p.s 
                          , nsim 
                          , plotsim){
  
  
  # Parameters for P_R (mdri.shape and mdri.scale), and unadjusted MDRI
  
  mdri.scale <- mdri.unadj/gamma(1+1/mdri.shape)
  
  fun.probrecent <- function(t, shape.in = mdri.shape, scale.in = mdri.scale){
    p <-  pweibull(t, shape = shape.in, scale = scale.in, lower.tail = FALSE)
    return(p)
  }

  mdri.unadj.rev <- integrate(fun.probrecent, lower = 0, upper = bigt)$value
  
  
  # Shape and scale parameters for Weibull distribution for time from testing to viral suppression
  # NB. Check limits of search space in function f.wblperc below
  
  vs.pars <- f.wblperc(vs.x1, vs.x2, vs.p1, vs.p2)
  vs.scale <- vs.pars$b 
  vs.shape <- vs.pars$k
  
  
  # Initialise plot
  
  par(mfrow = c(1,1))
  plot(0, 0, type = 'n', xlim = c(1, nsim+1), ylim = c(0, bigt), las = 1
       , ylab = 'Current MDRI approximation', xlab = 'Number of draws')
  abline(h = mdri.unadj.rev, col = 'red')
  abline(h = bigt, col = 'grey')
  grid()
  
  
  # Initial approximation of MDRI
  
  total.current.est <- 0 
  
  
  # For each draw of v...
  
  for (i in 1:nsim){
    
    
    # ...draw v_t2v...
    
    if (runif(1) <= p.vs.fixed) {
      v_t2v <- vs.fixed
    } else {
      v_t2v <- rweibull(1, shape = vs.shape, scale = vs.scale)
    }
    
    
    # ...draw v_i2t...
    
    grp <- (matrix(1:5, nrow =1) %*% rmultinom(1, 1, prob = c(p.i, p.r, p.f, p.n, p.s)))[1,1]
    if (grp == 1){ # i
      v_i2t <- rexp(1,rate = 1/delta.mean.i) + delaytotest.i
    } else if (grp == 2){ # r
      v_i2t <- runif(1, min = 0, max = delta.fixed.r)
    } else if (grp == 3){ # f
      delta.f <- runif(1, min = delta.min.f, max = delta.max.f)
      v_i2t <- runif(1, min = 0, max = delta.f)
    } else if (grp == 4){ # n
      v_i2t <- bigt + 10 
    } else { # s
      v_i2t <- 0
    }
    
    
    # ...integrate P_R until minimum of T and v...
    
    if (v_t2v + v_i2t >= bigt){
      add.new.est <- mdri.unadj.rev
    } else {
      add.new.est <- integrate(fun.probrecent, lower = 0, upper = v_t2v + v_i2t)$value 
    }
    
    
    # ...update approximation of MDRI...
    
    total.current.est <- total.current.est*(i-1)/i + 1/i*add.new.est
    
    
    # ...plot...
    
    if (i <= plotsim | i %% plotsim == 0 | i >= nsim - plotsim){
      # print(i)
      # print(add.new.est)
      # print(total.current.est)
      points(i,total.current.est, pch = 20, cex = 0.7)
    }
    
  }
  
  
  # Final MDRI value
  
  mdri.adj <-  total.current.est
  
  return(list(mdri.unadj.rev = mdri.unadj.rev
              , mdri.adj = mdri.adj
              , mdrired.perc = 100*(1-mdri.adj/mdri.unadj.rev)))
  
}



# Support functions -------------------------------------------------------


# Translate Weibull percentiles to shape and scale parameters

f.wblperc <- function(x1.fin, x2.fin, p1.fin, p2.fin
                      , mink.fin = 0, maxk.fin = 5 # limits of search space for shape parameter
                      , plot = FALSE){
  
  if (plot) {
    
    par(mfrow = c(1,2))
    
    plot(seq(mink.fin,maxk.fin,0.01)
         , f.tosolve.wblk(seq(mink.fin,maxk.fin,0.01)
                          , x1 = x1.fin, x2 = x2.fin, p1 = p1.fin, p2 = p2.fin)
         , type = 'l', xlab = 'Shape parameter', ylab = 'Function: to solve for root')
    abline(h=0, col = 'red')
    
  }
  
  wbl.k <- uniroot(f.tosolve.wblk,  interval = c(0, maxk.fin)
                   , x1 = x1.fin, x2 = x2.fin, p1 = p1.fin, p2 = p2.fin)$root # shape parameter
  wbl.b <- x1.fin/( log(1/(1-p1.fin)))^(1/wbl.k) # scale parameter 
  
  p1 <- 1-exp(-(x1.fin/wbl.b)^wbl.k) # Check probabilities
  p2 <- 1-exp(-(x2.fin/wbl.b)^wbl.k) # Check probabilities
  
  if (plot) {
    
    maxx <- qweibull(0.99, shape = wbl.k, scale = wbl.b)*1.2
    plot(seq(0, maxx, 0.01)
         , dweibull(seq(0, maxx, 0.01), shape = wbl.k, scale = wbl.b)
         , type = 'l', xlab = 'Time',  ylab = 'Weibull density function'
    )
    
    par(mfrow = c(1,1))
    
  }

    return(list(k = wbl.k, b = wbl.b, p1 = p1, p2 = p2))
  
}

f.tosolve.wblk <- function(k, x1, x2, p1, p2){
  y <- p2 - 1 + exp (-((x2/x1)^k*(log(1/(1-p1)))))
  return(y)
}



