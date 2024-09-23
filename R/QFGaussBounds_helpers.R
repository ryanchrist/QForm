####################
# Convenient function for directly calculating the target bounds with Quadrature (fails in the tails, hence the rest of these functions)
####################
RawQFGaussBounds <- function(cdf, f= "identity", max.abs.eta, sum.eta, sum.etasq, sum.eta.deltasq = 0, sum.etasq.deltasq = 0){

  nu <- 8*(sum.etasq.deltasq + (log(4)-1)*sum.etasq)
  L <- 1/(4*max.abs.eta)
  Er <- sum.eta.deltasq + sum.eta

  u <- function(y) ifelse(y <= 0, 1, ifelse(y < nu*L,exp(-0.5*(y^2)/nu),exp(0.5*nu*L^2-y*L) ))
  lower.integrand <- function(t,q){
    A <- cdf(t,density = T)
    ifelse(is.infinite(A),.Machine$double.xmax,A)*u(q - Er- t)}
  upper.integrand <- function(t,q){
    A <- cdf(t,density = T)
    ifelse(is.infinite(A),.Machine$double.xmax,A)*u(-(q - Er- t))}

  prelim.raw.bounds <- function(obs){
    lower <- max(0,cdf(obs-Er) - integrate(lower.integrand,lower = -Inf,upper=obs-Er,rel.tol = 1e-20,stop.on.error = FALSE,q = obs)$value)
    upper <- min(1,cdf(obs-Er) + integrate(upper.integrand,lower = obs-Er,upper=Inf,rel.tol = 1e-20,stop.on.error = FALSE, q = obs)$value)
    ans <- c(lower,upper,1-lower,1-upper)
    names(ans) <- c("lower.bound","upper.bound","one.minus.lower.bound","one.minus.upper.bound")
    ans
  }

  function(obs, parallel.sapply = sapply){
    t(parallel.sapply(obs,prelim.raw.bounds))
  }
  # Return function that takes in observed values
}


###################################
# Essential Integration Functions #
###################################

G1<-function(from,to,b){

  n1 <- length(from)
  n2 <- length(to)
  if(n1>1 & n2>1){stop("either from or to may be vectors but not both at the same time")}
  n <- max(n1,n2)

  # This performs the integral $ \int_\from^\to exp(-b*z) dz $ for any sign of b.
  #if(any(to<from)){warning("some to > from, so for those values returning -Inf")}
  if(any(is.infinite(from) & is.infinite(to))){stop("both from and psi cannot be infinite")}
  if(any(is.infinite(from) & from > 0)){stop("from cannot be positive infinity")}
  if(any(is.infinite(to) & to < 0)){stop("from cannot be negative infinity")}
  if(any(is.infinite(from)) & b > 0){stop("b cannot be positive when from is negative infinity --> integral does not converge")}
  if(any(is.infinite(to)) & b < 0){stop("b cannot be negative when to is infinity --> integral does not converge")}

  ifelse(to<=from, -Inf, ifelse(rep(b,n)==0,log(to-from),-log(abs(b))+log(-expm1(-abs(b)*(to-from)))-b*ifelse(rep(b,n)>0,from,to)))
}
#G1 is correct
#print(-G1(-20:-10,Inf,2)/log(10),digits=22)
#print(-log(integrate(f=function(z,b){exp(-b*z)},lower =-20,upper = Inf,b=2)$value)/log(10),digits=22)


G2 <- function(from, to, b, L){
  # This function performs Gauss Quadrature if L < 0.
  # either from or to may be vectors but not both simultaneously
  #if(any(to<from)){warning("some to > from, so for those values returning -Inf")}
  if(b>=1){stop("b must be < 1")}

  n1 <- length(to)
  n2 <- length(from)
  n <- max(n1,n2)
  ifelse(to<=from,-Inf,ifelse(rep(L,n)<0,
                              Vectorize(function(from,to,b,L){log(integrate(f=function(z,b,L){z^(-b)*exp(-L*z)},lower = from,upper = to,b=b,L=L,rel.tol = 1e-20,stop.on.error = FALSE)$value)},vectorize.args = c("from","to"))(from, to, b, L)
                              ,-(1-b)*log(L)+pgamma(L*to,1-b,log.p = TRUE)+lgamma(1-b) + (from >0)*log(-expm1(pgamma(L*from,1-b,log.p = TRUE)- pgamma(L*to,1-b,log.p = TRUE)))))
}
#G2 is correct
# print(G2(0,c(0.1,0.2),-3,-10),digits=22)
# print(log(integrate(f=function(z,b,L){z^(-b)*exp(-L*z)},lower = 0,upper = 0.2,b=-3,L=-10)$value),digits=22)



G3 <- function(from, to, mu, nu ){
  if(nu<=0){stop("nu must be > 0")}
  #if(any(to<from)){warning("some to > from, so for those values returning -Inf")}
  s<-sqrt(nu)
  from.tilde <- (from - mu)/s
  to.tilde <- (to - mu)/s
  ifelse(to<=from,-Inf,0.5*log(2*pi*nu) + pnorm(to.tilde,log.p = TRUE) + log(-expm1(pnorm(from.tilde,log.p = TRUE)-pnorm(to.tilde,log.p = TRUE))))
}
# #G3 is correct
#print(exp(G3(lower.prime,upper.prime,q-c2-b*nu,nu)),digits=22)
#print(integrate(f=function(z,mu,nu){exp(-0.5*(z-mu)^2/nu)},lower = lower.prime,upper = upper.prime,mu=q-c2-b*nu,nu=nu,rel.tol = 1e-13)$value,digits=22)
# print(exp(G3(-20,-10,0,.1)),digits=22)
# print(integrate(f=function(z,mu,nu){exp(-0.5*(z-mu)^2/nu)},lower = -20,upper = -10,mu=0,nu=.1,rel.tol = 1e-13)$value,digits=22)
#


G4 <- function(from, to, mu, nu ){
  if(nu<=0){stop("nu must be > 0")}
  #if(any(to<from)){warning("some to > from, so for those values returning -Inf")}
  s<-sqrt(nu)
  from.tilde <- (from - mu)/s
  to.tilde <- (to - mu)/s
  ifelse(to<=from,-Inf,0.5*log(2*pi*nu) + log( s*(dnorm(from.tilde)-dnorm(to.tilde))
                                               + mu*(pnorm(to.tilde)-pnorm(from.tilde))
  ))
}
# #G4 is correct
#print(exp(G4(lower.prime,upper.prime,q-c2-b*nu,nu)),digits=22)
#print(integrate(f=function(z,mu,nu){z*exp(-0.5*(z-mu)^2/nu)},lower = lower.prime,upper = upper.prime,mu=q-c2-b*nu,nu=nu,rel.tol = 1e-13)$value,digits=22)

# print(exp(G4(20,30,0,1)),digits=22)
# print(integrate(f=function(z,mu,nu){z*exp(-0.5*(z-mu)^2/nu)},lower = 20,upper = 30,mu=0,nu=1,rel.tol = 1e-13)$value,digits=22)

expG4 <- function(from, to, mu, nu ){
  if(nu<=0){stop("nu must be > 0")}
  #if(any(to<from)){warning("some to > from, so for those values returning -Inf")}
  s<-sqrt(nu)
  from.tilde <- (from - mu)/s
  to.tilde <- (to - mu)/s
  ifelse(to<=from,0,sqrt(2*pi*nu)*( s*(dnorm(from.tilde)-dnorm(to.tilde))
                                    + mu*(pnorm(to.tilde)-pnorm(from.tilde))
  ))
}


# G5 <- function(from, to, b, nu ){
#   if(any(from<0)){stop("from must be >= 0")}
#   if(any(to<from)){warning("some to > from, so for those values returning -Inf")}
#   if(b>=2){stop("b must be < 2")}
#   ifelse(to<=from,-Inf,-log(2) + G2(from^2,to^2,b/2,1/(2*nu)))
# }
# #G5 is correct
# print(exp(G5(-5,0:10,-2,2)),digits=22)
# print(integrate(f=function(z,b,nu){z^(1-b)*exp(-0.5*z^2/nu)},lower = 1,upper = 4,b=-20,nu=.02)$value,digits=22)



###################################
# WrapConcIneq for identity f     #
###################################

WrapConcIneq.identity <- function(c1,c2,nu,L){

  if(!all(is.numeric(c(c1,c2,nu,L)))){
    stop("c1,c2, nu, and L must all be numeric")
  }
  if(nu<=0){
    stop("nu must be positive")
  }
  if(L<=0){
    stop("L must be positive")
  }


  u = function(x) ifelse(x<0,1,ifelse(x < nu*L, exp(-0.5*(x^2)/nu), exp(0.5*nu*L^2-L*x)))
  l = function(x) ifelse(x>0,1,ifelse(x > -nu*L, exp(-0.5*(x^2)/nu), exp(0.5*nu*L^2+L*x)))

  list(
    "L" = L,
    "nu" = nu,
    "u" = u,
    "l" = l,
    "c1" = c1,
    "h1" = function(q,t){
      ifelse(t >= q-c1, 0 , ifelse(t > q-c1-nu*L, (q-c1-t)*exp(-0.5*((q-c1-t)^2)/nu)/nu,L*exp(0.5*nu*L^2-(q-c1)*L+L*t)))
    },
    "int.h1.const" = function(lower, upper, q, one.minus = FALSE){
      if(one.minus){
        if(lower>=upper){return(1)}
        res <- 1 - u(q-c1-upper) + u(q-c1-lower)
      }else{
        if(lower>=upper){return(0)}
        res <- u(q-c1-upper) - u(q-c1-lower)
      }
      res
    },
    "int.h1.expx" = function(lower, upper, q, a, b){
      if(lower>=upper){return(0)}
      lower.prime <- pmax(lower,q-c1-nu*L)
      upper.prime <- pmin(upper,q-c1)
      upper.prime.prime <- pmin(upper,q-c1-nu*L)

      A <- ifelse(lower.prime < upper.prime, max(min(exp((nu*b^2)/2 - b*(q-c1)-a),.Machine$double.xmax),.Machine$double.xmin)*((q-c1)*exp(G3(lower.prime,upper.prime,q-c1-b*nu,nu))-expG4(lower.prime,upper.prime,q-c1-b*nu,nu))/nu,0)
      B <- ifelse(lower < upper.prime.prime, L*exp(0.5*nu*L^2-(q-c1)*L-a+G1(lower,upper.prime.prime,b-L)), 0)
      A + B
    },
    "int.h1.explogx" = function(lower, upper, q, a, b){
      if(lower>=upper){return(0)}
      lower.prime <- pmax(lower,q-c1-nu*L)
      upper.prime <- pmin(upper,q-c1)
      upper.prime.prime <- pmin(upper,q-c1-nu*L)

      t.func <- function(t,q,c1,b,nu) { (q-c1-t)*exp(-0.5*((t-(q-c1))^2)/nu - b*log(t)-a-log(nu)) }

      A <- ifelse(lower.prime < upper.prime, Vectorize(function(f,lower,upper,q,c1,b,nu){integrate(f = f, lower = lower, upper = upper, q = q, c1 = c1, b = b, nu = nu,rel.tol = 1e-20,stop.on.error = F)$value},vectorize.args = c("lower","upper"))(t.func, lower.prime, upper.prime, q, c1, b, nu), 0)
      B <- ifelse(lower < upper.prime.prime, L*exp(0.5*nu*L^2-(q-c1)*L-a+G2(lower, upper.prime.prime, b, -L)), 0)
      A + B
    },
    "int.h1.explognegx" = function(lower, upper, q, a, b){
      if(lower>=upper){return(0)}
      lower.prime <- pmax(lower,q-c1-nu*L)
      upper.prime <- pmin(upper,q-c1)
      upper.prime.prime <- pmin(upper,q-c1-nu*L)

      z.func <- function(z,q,c1,b,nu) { (q-c1+z)*exp(-0.5*((z-(q-c1))^2)/nu - b*log(z) - a - log(nu))}

      A <- ifelse(lower.prime < upper.prime, Vectorize(function(f,lower,upper,q,c1,b,nu){integrate(f = f, lower = lower, upper = upper, q = q, c1 = c1, b = b, nu = nu,rel.tol = 1e-20,stop.on.error = F)$value},vectorize.args = c("lower","upper"))(z.func, -upper.prime, -lower.prime, q, c1, b, nu), 0)
      B <- ifelse(lower < upper.prime.prime, L*exp(0.5*nu*L^2-(q-c1)*L-a+G2(-upper.prime.prime,-lower, b, L)), 0)
      A + B
    },
    "c2" = c2,
    "h2" = function(q,t){
      ifelse(t <= q-c2, 0 , ifelse(t < q-c2+nu*L, (c2-q+t)*exp(-0.5*((c2-q+t)^2)/nu)/nu,L*exp(0.5*nu*L^2-(c2-q)*L-L*t)))
    },
    "int.h2.const" = function(lower, upper, q, one.minus = FALSE){
      if(one.minus){
        if(lower>=upper){return(1)}
        res <- 1 - u(c2-q+lower) + u(c2-q+upper)
      }else{
        if(lower>=upper){return(0)}
        res <- u(c2-q+lower) - u(c2-q+upper)
      }
      res
    },
    "int.h2.expx" = function(lower, upper, q, a, b){
      if(lower>=upper){return(0)}
      lower.prime <- pmax(lower,q-c2)
      upper.prime <- pmin(upper,q-c2+nu*L)
      lower.prime.prime <- pmax(lower,q-c2+nu*L)

      A <- ifelse(lower.prime < upper.prime, max(min(exp((nu*b^2)/2 - b*(q-c2)-a),.Machine$double.xmax),.Machine$double.xmin)*((c2-q)*exp(G3(lower.prime,upper.prime,q-c2-b*nu,nu))+expG4(lower.prime,upper.prime,q-c2-b*nu,nu))/nu, 0)
      B <- ifelse(lower.prime.prime < upper, L*exp(0.5*nu*L^2+(q-c2)*L-a+G1(lower.prime.prime,upper,b+L)), 0)
      A + B
    },
    "int.h2.explogx" = function(lower, upper, q, a, b){
      if(lower>=upper){return(0)}
      lower.prime <- pmax(lower,q-c2)
      upper.prime <- pmin(upper,q-c2+nu*L)
      lower.prime.prime <- pmax(lower,q-c2+nu*L)

      t.func <- function(t,q,c2,b,nu) { (c2-q+t)*exp(-0.5*((t-(q-c2))^2)/nu - b*log(t)-a-log(nu)) }

      A <- ifelse(lower.prime < upper.prime, Vectorize(function(f,lower,upper,q,c2,b,nu){integrate(f = f, lower = lower, upper = upper, q = q, c2 = c2, b = b, nu = nu,rel.tol = 1e-20,stop.on.error = F)$value},vectorize.args = c("lower","upper"))(t.func, lower.prime, upper.prime, q, c2, b, nu), 0)
      B <- ifelse(lower.prime.prime < upper, L*exp(0.5*nu*L^2+(q-c2)*L-a+G2(lower.prime.prime, upper, b, L)), 0)
      A + B
    },
    "int.h2.explognegx" = function(lower, upper, q, a, b){

      if(lower>=upper){return(0)}
      lower.prime <- pmax(lower,q-c2)
      upper.prime <- pmin(upper,q-c2+nu*L)
      lower.prime.prime <- pmax(lower,q-c2+nu*L)

      z.func <- function(z,q,c2,b,nu) { (c2-q-z)*exp(-0.5*((z-(c2-q))^2)/nu - b*log(z)-a-log(nu))}

      A <- ifelse(lower.prime < upper.prime, Vectorize(function(f,lower,upper,q,c2,b,nu){integrate(f = f, lower = lower, upper = upper, q = q, c2 = c2, b = b, nu = nu,rel.tol = 1e-20,stop.on.error = F)$value},vectorize.args = c("lower","upper"))(z.func, -upper.prime, -lower.prime, q, c2, b, nu), 0)
      B <- ifelse(lower.prime.prime < upper, L*exp(0.5*nu*L^2-(c2-q)*L-a+G2(-upper,-lower.prime.prime, b, -L)), 0)
      A + B
    }
  )
}


######################################
# Gauss Quadrature for "body" of CDF #
######################################

GaussQuadCDF<-function(from, to, lower.bound = TRUE, cdf,ql,qu,conc.ineqs){

  # This function integrates the "body" of the CDF against a function h.
  # Since h has a given mode that we don't want the Gauss quadrature routine to miss, we divide the integral into two regions of
  # integration, one leading up to and one leading away from the mode of h so that the the mode of h is always an endpoint of
  # the integration region passed to integrate (ensuring that the mode is not missed)

  components <- rep(0,2)

  c1 <- conc.ineqs$c1
  c2 <- conc.ineqs$c2
  nu <- conc.ineqs$nu

  if(lower.bound){
    if(from < (ql-c1-sqrt(nu)))
      components[1] <- integrate( function(t,ql){conc.ineqs$h1(ql,t)*cdf(t)}, lower=from,upper=min(to,ql-c1-sqrt(nu)),ql=ql,rel.tol = 1e-20,stop.on.error = FALSE)$value
    if(to > (ql-c1-sqrt(nu)) & from < (ql-c1) ){
      components[2] <- integrate( function(t,ql){conc.ineqs$h1(ql,t)*cdf(t)}, lower=max(ql-c1-sqrt(nu),from),upper=min(to,ql-c1),ql=ql,rel.tol=1e-20,stop.on.error = FALSE)$value
    }
  }

  if(!lower.bound){

    if((from < (qu-c2+sqrt(nu))) & (to > qu-c2) ){
      components[1] <- integrate( function(t,qu){conc.ineqs$h2(qu,t)*cdf(t)}, lower=max(from,qu-c2),upper=min(to,qu-c2+sqrt(nu)),qu=qu,rel.tol = 1e-20,stop.on.error = FALSE)$value
    }
    if(to > (qu-c2+sqrt(nu)))
      components[2] <- integrate( function(t,qu){conc.ineqs$h2(qu,t)*cdf(t)}, lower=max(from,qu-c2+sqrt(nu)),upper=to,qu=qu,rel.tol = 1e-20,stop.on.error = FALSE)$value
  }

  sum(components)
}


QFRestrictInterval <- function(x){
  x <- ifelse(x < 0, 0, x)
  ifelse(x > 1, 1, x)
}

QFSaferSum <- function(x){
  sum(x[order(abs(x),decreasing = F)])
}

