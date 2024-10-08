#' Bounds on the CDF of a Quadratic Form in Gaussians
#'
#' Returns a function for calculating upper and lower bounds on the CDF for random variables \eqn{Q_f = T_f + R_f} where only the CDF of \eqn{T_f} is known. These random variables have the form
#'
#' \deqn{T_f = \sum\limits_{i \in \mathcal{T}} f\left(\eta_i \right) A_i + \sigma Z_0}{T_f = \Sigma_{i \in T} f (\eta_i) A_i + \sigma Z_0},
#' \deqn{R_f = \sum\limits_{i \in \mathcal{R}} f\left(\eta_i \right) A_i}{T_f = \Sigma_{i \in R} f (\eta_i) A_i},
#' where each \eqn{A_i \sim \chi^2_{a_i}\left(\delta^2_i\right)}{A_i \sim \chi^2_{a_i}\left(\delta^2_i\right)} and \eqn{Z_0 \sim N(0,1)}{Z_0 \sim N(0,1)}, all mututally independent,
#' and \eqn{a_i = 1}{a_i = 1} for all \eqn{i \in \mathcal{R}}{i \in \mathcal{R}}. We aim to remove this final restriction in future work.
#'
#' If \code{max.abs.eta} < \code{.Machine$double.eps}, then the contribution of \eqn{R_f} to \eqn{Q_f} is ignored for numerical stability and the function returned is simply wrapper for the provided CDF of \eqn{T_f}.  If this is not desired, a user may want to consider rescaling \eqn{Q_f} to avoid this behavior.
#' Currently only \eqn{f = "identity"} is supported, but future versions will allow one to select \eqn{f} from a list or specify their own function with its corresponding bounds through a QFormFunction object.
#'
#' The returned bounds function takes a vector of observed values \code{q} at which to calculate bounds as it's main argument. If \code{q} is not known exactly, but only a lower bound \code{ql} and an upper bound \code{qu} are known, then those may provided instead of \code{q} and the returned bounds on the CDF will be valid for a \code{q} in \code{[ql,qu]}.  If \code{q} is provided, \code{ql} and \code{qu} are ignored.
#' The returned bounds function itself returns a matrix with four columns: \code{c("lower.bound","upper.bound","one.minus.lower.bound","one.minus.upper.bound")}.  The first and second columns are lower and upper bounds on the CDF at \code{q} respectively; the third and fourth columns are equal to one minus the first two columns but calculated separately by the function internally in order to maintain numerical precision in the upper tail.  Thus, it is strongly recommneded that users interested in upper tail p-values use the third and fourth columns rather than the first and second.
#'
#' The returned bounds function can also take a parallel version of sapply from a given R parallelization package via the optional argument \code{parallel.sapply}.  This can substantially speed up computation for long \code{q}.  See Example below and \code{\link{QFGauss}} for more details.
#'
#' \code{QFGaussBounds} cannot currently calculate bounds for a cdf returned by \code{QFGauss} that has a missing tail.  See \code{\link{QFGauss}} for more details.
#'
#' @param cdf function; the cdf of \eqn{T_f} returned by \code{QForm::QForm}
#' @param f character or QFormFunction object; the function \eqn{f} for the \eqn{Q_f} of interest.
#' @param max.abs.eta vector; element-wise upper bound on the absolute value of the \eqn{\eta_i} in \eqn{R_f} (see Details)
#' @param sum.eta vector; element-wise sum of the \eqn{\eta_i} in \eqn{R_f} (see Details)
#' @param sum.etasq vector; element-wise sum of the \eqn{\eta^2_i} in \eqn{R_f} (see Details)
#' @param sum.eta.deltasq vector; element-wise sum of the \eqn{\eta_i \delta^2_i} in \eqn{R_f} (see Details)
#' @param sum.etasq.deltasq vector; element-wise sum of the \eqn{\eta^2_i \delta^2_i} in \eqn{R_f} (see Details)
#' @param include.saddlepoint logical; if TRUE also return saddlepoint approximation based estimate of \eqn{Q_f} alongside bounds.  Currently only available when \eqn{f} = "identity." Default is FALSE.
#' @seealso \code{\link{QFGauss}}, \code{\link{TestQFGaussBounds}}
#' @return A vectorized function which evaluates upper and lower bounds on the CDF of \eqn{Q_f}.
#'
#' @examples
#'
#' f.eta <- c(-12, -7, 5, 7, -9, 10, 8)
#' delta2 <- c(2, 10, -4, 3, 8, -5, -12)^2
#'
#' cdf <- QFGauss(f.eta, delta2)
#'
#' bounds <- QFGaussBounds(cdf = cdf, f = "identity",
#'                         max.abs.eta = 10, sum.eta = 5, sum.etasq = 200)
#'
#' \dontrun{
#' # Evaluate the bounds at a set of points
#' xx <- seq(-1e3, 1e3, len = 6)
#' ## This may take 5 - 10 secs.
#' system.time(y <- bounds(xx))
#'
#' x <- seq(-1e3, 1e3, len = 1e3)
#' plot(x, cdf(x), type = "l", ylab = expression(CDF), xlab = expression(T[f]),
#'      main = expression(Bounds~on~CDF~of~T[f])) # CDF
#' points(xx, y[,1], col = "blue")
#' points(xx, y[,2], col = "red")
#'
#' # Generate diagnostic plots for bounds (currently TestQFGaussBounds only
#' # works for cases where the QFGauss produeced CDF has all df = 1.)
#' TestQFGaussBounds(QFGauss(c(1,5,-4,-3,10),c(2,-1,4,-5,5)^2),2)
#'
#' # The function returned by QFGaussBounds can be accelerated by passing it a
#' # parallel version of sapply.
#' # In this example we use only 2 parallel workers but more may be added
#' require(future.apply); plan(tweak(multiprocess,workers = 2))
#' system.time(y <- bounds(xx, parallel.sapply = future_sapply))
#' }
#'
#' @export
QFGaussBounds <- function(cdf, f = "identity", max.abs.eta, sum.eta, sum.etasq, sum.eta.deltasq = 0, sum.etasq.deltasq = 0, include.saddlepoint = FALSE){

  # This function takes a cdf generated by QFcdf and a function f with its associated arguments and returns
  # a function that returns bounds


  # Import tail features from the tail.features attribute of cdf
  tf <- attr(cdf,"tail.features")
  support <- tf$support
  ep.l <- tf$extrapolation.point.l
  ep.r <- tf$extrapolation.point.r
  a.l <- tf$a.l
  b.l <- tf$b.l
  a.r <- tf$a.r
  b.r <- tf$b.r


  # Don't compute bounds if one tail is missing (this may be relaxed in the future)
  if(any(is.na(c(a.l,b.l,a.r,b.r)))){stop("cdf cannot be bounded because at least one tail is missing: tail extrapolation in QFGauss failed, see ?QFGauss for details.")}

  if(max.abs.eta <= 0){stop("max.abs.eta must be positive")}

  if(max.abs.eta > sqrt(sum.etasq)){stop("max.abs.eta cannot be greater than sqrt(sum.etasq).")}

  if(max.abs.eta < .Machine$double.eps){
    warning("max.abs.eta is smaller than .Machine$double.eps so remainder term ignored (see Details).")

    bound.func <- function(X){
      ql <- X[1]
      qu <- X[2]

      ans <- c(cdf(ql), cdf(qu), cdf(ql, lower.tail = FALSE), cdf(qu, lower.tail = FALSE))

      names(ans) <- c("lower.bound","upper.bound","one.minus.lower.bound","one.minus.upper.bound")
      return(ans)
    }

  } else {

    # Create set of integration functions needed for integrating over the extrapolated tails of cdf corresponding to
    # the concentration bound for f
    if(f == "identity"){

      nu <- 8*(sum.etasq.deltasq + (log(4)-1)*sum.etasq)

      conc.ineqs <- WrapConcIneq.identity(sum.eta.deltasq + sum.eta,
                                          sum.eta.deltasq + sum.eta,
                                          nu,
                                          1/(4*max.abs.eta))
      saddlepoint.coeff <- 1/sqrt(2*pi*nu)

    }else{
      stop("For now, identity is the only valid option for f")
    }


    # Create bounding function
    bound.func <- function(X){
      ql <- X[1]
      qu <- X[2]

      # Initialize
      upper.components <- lower.components <- rep(0,5)
      # The first component here is a relic from a previous version.  That element will be removed in a later version.

      upper.components.4.alt <- lower.components.4.alt <- 1

      ### Compute Upper Bound

      if(qu-conc.ineqs$c2 < ep.l){

        if(support=="all.reals"){
          upper.components[2] <- conc.ineqs$int.h2.expx(qu-conc.ineqs$c2,ep.l,qu,a.l,b.l)
          upper.components[4] <- conc.ineqs$int.h2.const(ep.r,Inf,qu)
          upper.components.4.alt <- conc.ineqs$int.h2.const(ep.r,Inf,qu, one.minus = TRUE)
          upper.components[5] <- -conc.ineqs$int.h2.expx(ep.r,Inf,qu,a.r,b.r)
        }
        if(support=="pos.reals"){
          upper.components[2] <- conc.ineqs$int.h2.explogx(max(qu-conc.ineqs$c2,0),ep.l,qu,a.l,b.l)
          upper.components[4] <- conc.ineqs$int.h2.const(ep.r,Inf,qu)
          upper.components.4.alt <- conc.ineqs$int.h2.const(ep.r,Inf,qu, one.minus = TRUE)
          upper.components[5] <- -conc.ineqs$int.h2.expx(ep.r,Inf,qu,a.r,b.r)
        }
        if(support=="neg.reals"){
          upper.components[2] <- conc.ineqs$int.h2.expx(qu-conc.ineqs$c2,ep.l,qu,a.l,b.l)
          upper.components[4] <- conc.ineqs$int.h2.const(ep.r,Inf,qu) # This integrates out to infinity because we need to account for the constant plateau of the NSD density above 0
          upper.components.4.alt <- conc.ineqs$int.h2.const(ep.r,Inf,qu, one.minus = TRUE)
          upper.components[5] <- -conc.ineqs$int.h2.explognegx(ep.r,0,qu,a.r,b.r)
        }

        #upper.components[3] <- boole(ep.l,t.cdf$x, conc.ineqs$h2(qu,t.cdf$x)*t.cdf$y) #Boole integral from ep.l to t.cdf$x[n]
        upper.components[3] <- GaussQuadCDF(ep.l, ep.r, F, cdf, ql, qu, conc.ineqs) #Gauss Quad integral from ep.l to t.cdf$x[n]

      }

      if((qu-conc.ineqs$c2 >= ep.l) & (qu-conc.ineqs$c2 < ep.r)){

        if(support=="all.reals"){
          upper.components[4] <- conc.ineqs$int.h2.const(ep.r,Inf,qu)
          upper.components.4.alt <- conc.ineqs$int.h2.const(ep.r,Inf,qu, one.minus = TRUE)
          upper.components[5] <- -conc.ineqs$int.h2.expx(ep.r,Inf,qu,a.r,b.r)
        }
        if(support=="pos.reals"){
          upper.components[4] <- conc.ineqs$int.h2.const(ep.r,Inf,qu)
          upper.components.4.alt <- conc.ineqs$int.h2.const(ep.r,Inf,qu, one.minus = TRUE)
          upper.components[5] <- -conc.ineqs$int.h2.expx(ep.r,Inf,qu,a.r,b.r)
        }
        if(support=="neg.reals"){
          upper.components[4] <- conc.ineqs$int.h2.const(ep.r,Inf,qu) # This integrates out to infinity because we need to account for the constant plateau of the NSD density above 0
          upper.components.4.alt <- conc.ineqs$int.h2.const(ep.r,Inf,qu, one.minus = TRUE)
          upper.components[5] <- -conc.ineqs$int.h2.explognegx(ep.r,0,qu,a.r,b.r)
        }

        #upper.components[3] <- boole(qu,t.cdf$x, conc.ineqs$h2(qu,t.cdf$x)*t.cdf$y) #Boole integral from qu to t.cdf$x[n]
        upper.components[3] <- GaussQuadCDF(qu-conc.ineqs$c2, ep.r, F, cdf, ql, qu, conc.ineqs) #Gauss Quad integral from qu to t.cdf$x[n]

      }

      if(qu-conc.ineqs$c2 >= ep.r){

        if(support=="all.reals"){
          upper.components[4] <- conc.ineqs$int.h2.const(qu-conc.ineqs$c2,Inf,qu)
          upper.components.4.alt <- conc.ineqs$int.h2.const(qu-conc.ineqs$c2,Inf,qu, one.minus = TRUE)
          upper.components[5] <- -conc.ineqs$int.h2.expx(qu-conc.ineqs$c2,Inf,qu,a.r,b.r)
        }
        if(support=="pos.reals"){
          upper.components[4] <- conc.ineqs$int.h2.const(qu-conc.ineqs$c2,Inf,qu)
          upper.components.4.alt <- conc.ineqs$int.h2.const(qu-conc.ineqs$c2,Inf,qu, one.minus = TRUE)
          upper.components[5] <- -conc.ineqs$int.h2.expx(qu-conc.ineqs$c2,Inf,qu,a.r,b.r)
        }
        if(support=="neg.reals"){
          upper.components[4] <- conc.ineqs$int.h2.const(qu-conc.ineqs$c2,Inf,qu) # This integrates out to infinity because we need to account for the constant plateau of the NSD density above 0
          upper.components.4.alt <- conc.ineqs$int.h2.const(qu-conc.ineqs$c2,Inf,qu, one.minus = TRUE)
          upper.components[5] <- -conc.ineqs$int.h2.explognegx(qu-conc.ineqs$c2,0,qu,a.r,b.r)
        }

      }


      ### Compute Lower Bound

      if(ql-conc.ineqs$c1 > ep.r){

        if(support=="all.reals"){
          lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,ep.l,ql,a.l,b.l)
          lower.components[4] <- conc.ineqs$int.h1.const(ep.r,ql-conc.ineqs$c1,ql)
          lower.components.4.alt <- conc.ineqs$int.h1.const(ep.r,ql-conc.ineqs$c1,ql, one.minus = TRUE)
          lower.components[5] <- -conc.ineqs$int.h1.expx(ep.r,ql-conc.ineqs$c1,ql,a.r,b.r)
        }
        if(support=="pos.reals"){
          lower.components[2] <- conc.ineqs$int.h1.explogx(0,ep.l,ql,a.l,b.l)
          lower.components[4] <- conc.ineqs$int.h1.const(ep.r,ql-conc.ineqs$c1,ql)
          lower.components.4.alt <- conc.ineqs$int.h1.const(ep.r,ql-conc.ineqs$c1,ql, one.minus = TRUE)
          lower.components[5] <- -conc.ineqs$int.h1.expx(ep.r,ql-conc.ineqs$c1,ql,a.r,b.r)
        }
        if(support=="neg.reals"){
          lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,ep.l,ql,a.l,b.l)
          lower.components[4] <- conc.ineqs$int.h1.const(ep.r,ql-conc.ineqs$c1,ql) # This integrates out to ql instead of min(ql,0) because we need to account for the constant plateau of the NSD density above 0
          lower.components.4.alt <- conc.ineqs$int.h1.const(ep.r,ql-conc.ineqs$c1,ql, one.minus = TRUE)
          lower.components[5] <- -conc.ineqs$int.h1.explognegx(ep.r,min(0,ql-conc.ineqs$c1),ql,a.r,b.r)
        }

        #lower.components[3] <- boole(ep.l,t.cdf$x, conc.ineqs$h1(ql,t.cdf$x)*t.cdf$y) #Boole integral from ep.l to t.cdf$x[n]
        lower.components[3] <- GaussQuadCDF(ep.l, ep.r, TRUE, cdf, ql, qu, conc.ineqs) #Gauss Quad integral from ep.l to t.cdf$x[n]

      }

      if((ql-conc.ineqs$c1 > ep.l) & (ql-conc.ineqs$c1 <= ep.r)){

        if(support=="all.reals"){
          lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,ep.l,ql,a.l,b.l)
        }
        if(support=="pos.reals"){
          lower.components[2] <- conc.ineqs$int.h1.explogx(0,ep.l,ql,a.l,b.l)
        }
        if(support=="neg.reals"){
          lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,ep.l,ql,a.l,b.l)
        }

        #lower.components[3] <- boole(ql,t.cdf$x, conc.ineqs$h1(ql,t.cdf$x)*t.cdf$y,int.to.right = F) #Boole integral from ep.l to t.cdf$x[n]
        lower.components[3] <- GaussQuadCDF(ep.l, ql-conc.ineqs$c1, TRUE, cdf, ql, qu, conc.ineqs) #Gauss Quad integral from ep.l to ql
      }

      if(ql-conc.ineqs$c1 <= ep.l){

        if(support=="all.reals"){
          lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,ql-conc.ineqs$c1,ql,a.l,b.l)
        }
        if(support=="pos.reals"){
          lower.components[2] <- conc.ineqs$int.h1.explogx(0,ql-conc.ineqs$c1,ql,a.l,b.l)
        }
        if(support=="neg.reals"){
          lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,ql-conc.ineqs$c1,ql,a.l,b.l)
        }

      }

      # Sum across components

      one.minus.upper.components <- - upper.components
      one.minus.upper.components[4] <- upper.components.4.alt

      one.minus.lower.components <- - lower.components
      one.minus.lower.components[4] <- lower.components.4.alt

      #RC: Could revisit this approach to numerically adding in the two parts of the const integral, but that seems to be
      # the part that causes the greatest trouble and this seems to work well for now so will avoid introducing other dependencies
      # though using a more standard safe sum approach on the components plus the two coming from the const integral might be
      # better as we enter either tail of the distribution

      ans <- QFRestrictInterval(c(QFSaferSum(lower.components),
                                  QFSaferSum(upper.components),
                                  QFSaferSum(one.minus.lower.components),
                                  QFSaferSum(one.minus.upper.components)))

      if(include.saddlepoint & f == "identity"){
        ans <- c(ans[1], saddlepoint.coeff * (ans[1] + ans[2]),ans[2],
                 ans[3], saddlepoint.coeff * (ans[3] + ans[4]),ans[4])
        names(ans) <- c("lower.bound","sp.approx","upper.bound","one.minus.lower.bound","one.minus.sp.approx","one.minus.upper.bound")
        return(ans)
      }else{
        names(ans) <- c("lower.bound","upper.bound","one.minus.lower.bound","one.minus.upper.bound")
        return(ans)
      }
    }
  }


  wrapped.bound.func <- function(q = NULL, ql = NULL, qu = NULL, parallel.sapply = base::sapply){
    # If q is specified, it over-rides qu and ql.
    if(!is.null(q)){
      if(!is.numeric(q)){stop("q must be numeric")}
      return(t(parallel.sapply(split(cbind(q,q),1:length(q)),bound.func)))
    }else{
      if( !(is.numeric(ql) & is.numeric(qu))  ){stop("ql and qu must be numeric")}
      if(length(ql)!=length(qu)){stop("ql and qu must have the same length")}
      if(any(qu<ql)){stop("ql must be less than or equal to qu")}
      return(t(parallel.sapply(split(cbind(ql,qu),1:length(ql)),bound.func)))
    }
  }

  class(wrapped.bound.func) <- c("QFGaussBoundsFunc",class(wrapped.bound.func))
  wrapped.bound.func

}


#' Test function for a QFGaussBounds
#'
#' Compares the CDF bounds inferred by QFGaussBounds to a truncated approximation of the CDF and a naive quadrature-based implementation of the bounds.
#'
#' Here, \code{fullcdf} is taken to be the CDF of the target random variable \eqn{Q_f} (see documentation of \code{\link{QFGauss}} for definitions).
#' Four plots are produced.  The top-left plot overlays the target CDF of \eqn{Q_f}, \code{fullcdf}, computed by \code{QFGauss} (in black) and a truncated approximation to that CDF (in orange) based on simply adding the expectation of the remainder term \eqn{R_f} to the top-k truncated version of \code{fullcdf}, \eqn{T_f}.  By "top-k" here we mean taking the terms of \eqn{Q_f} with the largest magnitude coefficients, \code{f.eta}, and using that to define \eqn{T_f}, which is what is done in \code{TestQFGaussBounds} internally.
#' The green line is a similar approximation but where \eqn{R_f} is approximated with a moment-matching gaussian.
#' The upper and lower bounds on \code{fullcdf} computed by \code{QFGaussBounds} are plotted as red and blue circles respectively.
#' The upper and lower bounds on \code{fullcdf} computed by a naive quadrature-based implementation of the bounds are plotted as red and blue Xs respectively.
#' The top-right plot shows the difference between the truncated approximation of the CDF and \code{fullcdf} in log space.  It may be interpreted as follows.  The x-axis plots the -log_10 p-value one would have reported based on the truncated approximation alone.  The y-axis is the difference between
#' the true -log_10 p-value and the approximate -log_10 p-value.  The difference in p-values in the upper tail is plotted with a solid line.  The difference in p-values in the lower tail is plotted with a dashed line.  This plot effectively shows how far one might be misled by the truncated approximation shifted by the expectation of the remainder terms.
#' The two bottom plots allow comparison of the empirical CDF (in red) with the computed CDF (in black) in each tail.
#'
#' @param fullcdf QFGaussCDF; the target CDF including all terms; currently TestQFGaussBounds only works for cases where the QFGauss produeced CDF has all df = 1.
#' @param k numeric; the number of truncated terms provided to QFGaussBounds from which to bound fullcdf
#' @param n.bound.points numeric; the number of points at which to evaluate the bounds for plotting
#' @param lower.tail.end numeric; the -log_10 lower tail probability at which to start each x-axis (default = 20)
#' @param upper.tail.end numeric; the -log_10 upper tail probability at which to end each x-axis (default = 20 )
#' @param parallel.sapply function; a user-provided version of \code{sapply}, see Details.
#' @return There is nothing returned.
#' @seealso \code{\link{QFGaussBounds}}, \code{\link{TestQFGauss}}
#' @examples
#' TestQFGaussBounds(QFGauss(c(1,5,-4,-3),c(2,-1,4,-5)^2),2)
#'
#' @importFrom graphics par points lines abline
#' @export
TestQFGaussBounds <- function(fullcdf, k = min(20,floor(length(attr(fullcdf,"f.eta"))/2)), n.bound.points = 16,
                              lower.tail.end = 20, upper.tail.end = 20, parallel.sapply = base::sapply){

  message("Comparing the provided fullcdf to a truncated approximation base on the first k = ",k," terms.  This may take a minute.")
  if(class(fullcdf)[1]!="QFGaussCDF"){stop("fullcdf must be of class QFGaussCDF")}

  fft_used <- attr(fullcdf,"fft_used")
  f.eta <- attr(fullcdf,"f.eta")
  if(length(f.eta)==1){stop("fullcdf cannot be compared to truncated version because length(f.eta) is only 1.")}

  # Don't compute bounds if any of the df are not equal to 1.
  if(!all(attr(fullcdf,"df")==1)){stop("TestQFGaussBounds currently only supported for cases where all df == 1")}

  delta2 <- attr(fullcdf,"delta2")
  mu <- attr(fullcdf,"mu")
  Q.sd <- attr(fullcdf,"Q.sd")

  tf <- attr(fullcdf,"tail.features")
  support <- tf$support
  ep.l <- tf$extrapolation.point.l
  ep.r <- tf$extrapolation.point.r
  a.l <- tf$a.l
  b.l <- tf$b.l
  a.r <- tf$a.r
  b.r <- tf$b.r


  if(support == "all.reals"){
    x.max <-uniroot(function(z) {- fullcdf(z,lower.tail = F,log.p = T) / log(10) - upper.tail.end},
                    lower = ep.r, upper = ep.r + 0.1*Q.sd,tol = .Machine$double.eps, extendInt = "upX")$root
    x.min <-uniroot(function(z) {- fullcdf(z,lower.tail = T,log.p = T) / log(10) - lower.tail.end},
                    lower = ep.l - 0.1*Q.sd, upper = ep.l,tol = .Machine$double.eps, extendInt = "downX")$root
  }

  if(support == "pos.reals"){
    x.min <- 0
    x.max <-uniroot(function(z) {- fullcdf(z,lower.tail = F,log.p = T) / log(10) - upper.tail.end},
                    lower = ep.r, upper = ep.r + 0.1*Q.sd,tol = .Machine$double.eps, extendInt = "upX")$root
  }

  if(support == "neg.reals"){
    x.min <-uniroot(function(z) {- fullcdf(z,lower.tail = T,log.p = T) / log(10) - lower.tail.end},
                    lower = ep.l - 0.1*Q.sd, upper = ep.l,tol = .Machine$double.eps, extendInt = "downX")$root
    x.max <- 0
  }

  x <- seq(x.min,x.max,len=1e5)
  xx <- seq(x.min,x.max,len=n.bound.points)

  old.par <- par(no.readonly = T)
  par(mfrow=c(2,2))


  tcdf <- QFGauss(f.eta[1:k],delta2[1:k],parallel.sapply = parallel.sapply)

  max.abs.eta <- abs(f.eta[k])
  sum.eta <- sum(f.eta[(k+1):length(f.eta)])
  sum.etasq <- sum(f.eta[(k+1):length(f.eta)]^2)
  sum.eta.deltasq <- sum(f.eta[(k+1):length(f.eta)]*delta2[(k+1):length(f.eta)])
  sum.etasq.deltasq <- sum((f.eta[(k+1):length(f.eta)]^2)*delta2[(k+1):length(f.eta)])

  # if f.eta[k] is really large in magnitude compared to the f.eta[(k+1):length(f.eta)] , it's
  # possible that abs(f.eta[k]) does not place a good enough bound on max.abs.eta.  Here, we use the
  # fact that the l_2 norm must always dominate the l_\infty to obtain a better bound for max.abs.eta
  if(max.abs.eta > sqrt(sum.etasq)){max.abs.eta <- sqrt(sum.etasq)}

  bound.func <- QFGaussBounds(tcdf,
                              "identity",
                              max.abs.eta,
                              sum.eta,
                              sum.etasq,
                              sum.eta.deltasq,
                              sum.etasq.deltasq)

  raw.bound.func <- RawQFGaussBounds(tcdf,
                                     "identity",
                                     max.abs.eta,
                                     sum.eta,
                                     sum.etasq,
                                     sum.eta.deltasq,
                                     sum.etasq.deltasq)

  bounds <- bound.func(xx,parallel.sapply = parallel.sapply)
  raw.bounds <- raw.bound.func(xx,parallel.sapply = parallel.sapply)


  Er <- sum.eta.deltasq + sum.eta
  approxfullcdf <- function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) tcdf(x-Er, density = density, lower.tail = lower.tail, log.p = log.p)
  attr(approxfullcdf,"mu") <- attr(tcdf,"mu") + Er
  attr(approxfullcdf,"Q.sd") <- attr(tcdf,"Q.sd")


  VARr <- 2*sum.etasq + 4*sum.etasq.deltasq

  gauss.tcdf <- QFGauss(f.eta[1:k],delta2[1:k],sigma = sqrt(VARr), parallel.sapply = parallel.sapply)

  gauss.approxfullcdf <- function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) gauss.tcdf(x-Er, density = density, lower.tail = lower.tail, log.p = log.p)
  attr(gauss.approxfullcdf,"mu") <- attr(gauss.tcdf,"mu") + Er
  attr(gauss.approxfullcdf,"Q.sd") <- attr(gauss.tcdf,"Q.sd")


  old.par <- par(no.readonly = T)
  par(mfrow=c(2,2))

  plot(x,fullcdf(x),type="l",lwd=1.5,
       ylab = expression(CDF),xlab=expression(Q[f]), main=expression(CDF~of~Q[f])) # CDF as calculated by QForm
  points(xx,bounds[,1],col="blue")
  points(xx,bounds[,2],col="red")
  points(xx,raw.bounds[,1],col="blue",pch=4)
  points(xx,raw.bounds[,2],col="red",pch=4)
  lines(x,approxfullcdf(x),type="l",col="darkorange",lwd=1.5)
  lines(x,gauss.approxfullcdf(x),type="l",col="green",lwd=1.5)



  true.pval.uppertail <- fullcdf(suppressWarnings(uniroot(function(z) {- approxfullcdf(z,lower.tail = F,log.p = T) / log(10) - 8},
                                                          lower = attr(approxfullcdf,"mu")-2*attr(approxfullcdf,"Q.sd"),
                                                          upper = attr(approxfullcdf,"mu") + 0.1*attr(approxfullcdf,"Q.sd"),
                                                          tol = .Machine$double.eps, extendInt = "upX")$root),lower.tail = F)

  true.pval.lowertail <- fullcdf(suppressWarnings(uniroot(function(z) {- approxfullcdf(z,log.p = T) / log(10) - 8},
                                                          lower = attr(approxfullcdf,"mu")-0.1*attr(approxfullcdf,"Q.sd"),
                                                          upper = attr(approxfullcdf,"mu") + 2*attr(approxfullcdf,"Q.sd"),
                                                          tol = .Machine$double.eps, extendInt = "downX")$root))

  yy.uppertail <- -fullcdf(x,lower.tail = F, log.p = T)/log(10) + approxfullcdf(x,lower.tail = F, log.p = T)/log(10)
  yy.lowertail <- -fullcdf(x, log.p = T)/log(10) + approxfullcdf(x, log.p = T)/log(10)
  xx.uppertail <- -approxfullcdf(x,lower.tail = F, log.p = T)/log(10)
  xx.lowertail <- -approxfullcdf(x, log.p = T)/log(10)

  gauss.yy.uppertail <- -fullcdf(x,lower.tail = F, log.p = T)/log(10) + gauss.approxfullcdf(x,lower.tail = F, log.p = T)/log(10)
  gauss.yy.lowertail <- -fullcdf(x, log.p = T)/log(10) + gauss.approxfullcdf(x, log.p = T)/log(10)
  gauss.xx.uppertail <- -gauss.approxfullcdf(x,lower.tail = F, log.p = T)/log(10)
  gauss.xx.lowertail <- -gauss.approxfullcdf(x, log.p = T)/log(10)

  plot(xx.uppertail,
       yy.uppertail,
       type = "l", lwd = 1.5, ylab = "Approx. - True -log10 p-values",
       xlab = "Approx. -log10 p-value",
       ylim = range(pretty(c(yy.uppertail,yy.lowertail))),
       xlim = range(pretty(c(xx.uppertail,xx.lowertail))),
       main = paste("At Naive p = 1e-8, upper tail p =",
                    signif(true.pval.uppertail,digits=3),", lower tail p =",signif(true.pval.lowertail,digits=3)),
       font.main = 1, col="darkorange"
  )

  lines(xx.lowertail,
        yy.lowertail,
        type = "l", lwd = 1.5,lty=3, col="darkorange")

  lines(gauss.xx.uppertail,
        gauss.yy.uppertail,
        type = "l", lwd = 1.5, col="green")

  lines(gauss.xx.lowertail,
        gauss.yy.lowertail,
        type = "l", lwd = 1.5,lty=3, col="green")
  abline(h=0,lty=4)

  plot(x,-fullcdf(x,log.p = T)/log(10),type="l",lwd=1.5,
       ylab = expression(-log[10](CDF)),main=expression(-log[10](CDF)), xlab=expression(Q[f]))# CDF as calculated by QForm
  points(xx,-log10(bounds[,1]),col="blue")
  points(xx,-log10(bounds[,2]),col="red")
  points(xx,-log10(raw.bounds[,1]),col="blue",pch=4)
  points(xx,-log10(raw.bounds[,2]),col="red",pch=4)
  lines(x,-approxfullcdf(x,log.p = T)/log(10),type="l",col="darkorange",lwd=1.5)
  lines(x,-gauss.approxfullcdf(x,log.p = T)/log(10),type="l",col="green",lwd=1.5)


  plot(x,-fullcdf(x,lower.tail = F,log.p = T)/log(10),type="l",lwd=1.5,
       ylab = expression(-log[10](1 - CDF)),main=expression(-log[10](1-CDF)), xlab=expression(Q[f])) # CDF as calculated by QForm
  points(xx,-log10(bounds[,3]),col="blue")
  points(xx,-log10(bounds[,4]),col="red")
  points(xx,-log10(raw.bounds[,3]),col="blue",pch=4)
  points(xx,-log10(raw.bounds[,4]),col="red",pch=4)
  lines(x,-approxfullcdf(x,lower.tail = F,log.p = T)/log(10),type="l",col="darkorange",lwd=1.5)
  lines(x,-gauss.approxfullcdf(x,lower.tail = F,log.p = T)/log(10),type="l",col="green",lwd=1.5)


  par(old.par)
}

