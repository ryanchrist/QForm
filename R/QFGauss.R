#' Fast CDF/PDF of a Quadratic Form in Gaussians
#'
#' Returns the CDF and PDF for random variables \eqn{T_f}{T_f} of the form \deqn{T_f = \sum\limits_i f\left(\eta_i \right) \left(Z_i + \delta_i)^2 + \sigma Z_0}{T_f = \Sigma_i f (\eta_i) (Z_i + \delta_i)^2 + \sigma Z_0} where \eqn{Z_i \sim N(0,1).}{Z_i ~ N(0,1).}
#'
#' By using the fast Fourier transform and various adjustments for numerical precision, this function is faster and more reliable than Davie's method and related approaches, especially when the returned CDF or PDF is to be evaluated at many points.
#'
#' The returned function has three optional, logical arguments.  The first is a \code{density}, which when \code{TRUE}, prompts the function to evaluate the PDF rather than the CDF.  \code{density} defaults to FALSE.  \code{lower.tail} returns 1 minus the CDF when \code{TRUE} (not used if \code{density}==TRUE) and is highly recommended for those interested in the upper tail of \eqn{T_f}.  \code{log.p} returns the desired probabilities in log space.
#'
#' \code{parallel.sapply}, by default, is set to \code{base:sapply}.  However, it allows the user to supply a parallelized version of \code{sapply} (eg: \code{future_sapply} from the \code{future.apply} package) to help speed up the calculation of the CDF.  This is helpful in cases where \code{length(f.eta)} is large.
#'
#' \code{n} is the number of sub-intervals used in the left-sided Reimann integral approximation of the Fourier transform carried out by \code{stats::fft}.  The default 2^16-1 should work for the vast majority of cases, but n may need to be increased to achieve accurate CDF estimation when \eqn{T_f} has many terms (when \code{f.eta} is long).
#'
#' Since \code{stats::fft} can only evaluate the CDF up to double precision, we extrapolate the tails of \eqn{T_f}.  QForm automatically detects the region where the estimated CDF begins to lose precision.  A log-linear function is used for tails that go out to infinity and a function of the form \eqn{\alpha |x|^\beta} is used for tails truncated at 0 (when all of the \code{f.eta} have the same sign).  These extrapolated tails, motivated by the form of the characteristic function, provide accurate approximations in most cases when compared against a quad-precision implementation (not yet included in \code{QForm}).
#'
#' Our current tail extrapolation scheme can become unstable or fail in cases where the target distribution is extremely skewed.  In these cases, one of the tails decays too rapidly to be estimated with the given number of FFT grid points (set by QFGauss optional argument n).  \code{QFGaussBounds} cannot currently calculate bounds for a cdf returned by \code{QFGauss} that has a missing tail.  While we plan to address extremely skewed cases in future versions by deploying a second FFT when needed, for now, we recommend that users who really care about estimation of the thin tail
#' or obtaining bounds with \code{QFGaussBounds} to try increasing the number of FFT grid points, \code{n}, passed to \code{QFGauss}.
#'
#' A note on unbounded densities: The density of \eqn{T_f} is guaranteed to be bounded if length(f.eta) > 2 and there is no trouble in density estimation posed by asymptotes.  In the length(f.eta)==2 case, if the two components of f.eta are of opposite signs, then the density of \eqn{T_f}{T_f} may have an asymptote at some value \eqn{t}{t}.  While the density in the neighborhood around that \eqn{t} should be accurately calculated, due to the FFT and spline interpolation approach used, the density reported at \eqn{t}{t} may be reported as some finite rather than as Inf.  In the length(f.eta)==1 case, QFGauss resorts to dchisq and the density at 0 is accurately reported as Inf.
#'
#'
#' @param f.eta vector; real-valued coefficients, \eqn{f(\eta_i)}, (may be positive or negative)
#' @param delta2 vector; non-negative real-valued non-centrality parameters for each term (default is 0s). As is standard for chi-squared non-centrality parameters, these are assumed to be already summed across terms when df > 1.
#' @param df vector; positive real-valued degrees of freedom for each term (default is vector of 1s), can really increase speed if there are many redundant terms
#' @param sigma numeric; standard deviation of optional Gaussian term \eqn{Z_0}, default is 0 (no Gaussian term added)
#' @param n integer; number of points at which to evaluate the characteristic function of \eqn{T_f}, must be odd (see Details).
#' @param parallel.sapply function; a user-provided version of \code{sapply}, see Details.
#'
#' @return A function that evaluates the CDF or PDF of \eqn{T_f}.
#' @seealso \code{\link{QFGaussBounds}}, \code{\link{TestQFGauss}}
#' @examples
#' f.eta <- c(-12, -7, 5, 7, -9, 10, 8)
#' delta2 <- c(2, 10, -4, 3, 8, -5, -12)^2
#' df <- c(1.1,5.2,0.4,10,1,2.5,1)
#'
#' cdf <- QFGauss(f.eta, delta2, df)
#'
#' # Inspect computed CDF
#' plot(cdf)
#'
#' # Plot computed CDF at desired points
#' x <- seq(-1500, 2000, len = 1e3)
#' plot(x,cdf(x),type="l", ylab = expression(CDF),xlab=expression(T[f]), main=expression(CDF~of~T[f])) # CDF
#' plot(x,cdf(x,density = TRUE),type="l", ylab = expression(PDF),xlab=expression(T[f]), main=expression(PDF~of~T[f])) # PDF
#'
#' # Compare computed CDF to empirical CDF of target distribution based on 10,000 samples
#' TestQFGauss(cdf)
#'
#' # QFGauss can be accelerated by passing it a parallel version of sapply
#' \dontrun{
#' # In this example we use only 2 parallel workers but more may be added
#' require(future.apply); plan(tweak(multiprocess,workers = 2))
#' f.eta <- 5 * rnorm(500)
#' system.time(cdf <- QFGauss(f.eta))
#' system.time(cdf <- QFGauss(f.eta, parallel.sapply = future_sapply))
#' }
#'
#'
#'
#' @export
QFGauss <- function(f.eta, delta2 = rep(0,length(f.eta)),  df = rep(1,length(f.eta)), sigma = 0, n = 2^16-1, parallel.sapply = base::sapply){

  if(n%%2==0){stop("n must be odd")}

  if(length(f.eta)!=length(delta2)){stop("delta2 does not have the same length as f.eta.")}
  if(length(f.eta)!=length(df)){stop("df does not have the same length as f.eta.")}
  if(length(sigma)!=1){stop("sigma must have length 1")}

  if(!all(is.finite(f.eta))){stop("All f.eta must be finite.")}
  if(!all(is.finite(delta2))){stop("All delta2 must be finite.")}
  if(!all(is.finite(df))){stop("All df must be finite.")}
  if(!is.finite(sigma)){stop("sigma must be finite.")}

  if(all(f.eta==0)){stop("All f.eta are zero.")}
  if(any(delta2 < 0)){stop("All delta2 must be >= 0")}
  if(any(df <= 0)){stop("All df must be > 0")}
  if(sigma < 0){stop("sigma cannot be negative")}

  # Reorder f.eta so that the largest magnitude f.eta come first.  Permute delta2 and df accordingly.
  f.eta.order <- order(abs(f.eta),decreasing=T)
  f.eta <- f.eta[f.eta.order]
  delta2 <- delta2[f.eta.order]
  df <- df[f.eta.order]


  if(all(f.eta==f.eta[1]) & sigma==0){

    total_df <- sum(df)
    C <- abs(f.eta[1])
    total_ncp <- sum(delta2)

      cdf.func <- function(q, density = FALSE, lower.tail = TRUE, log.p = FALSE){
        if(density){
          if(log.p){
            return( dchisq((2*(f.eta[1]>0)-1)*q/C,total_df,total_ncp,TRUE)-log(C) )
          }else{
            return( dchisq((2*(f.eta[1]>0)-1)*q/C,total_df,total_ncp)/C )
          }
        }else{
          if(f.eta[1] > 0){
            return( pchisq(q/C,total_df,total_ncp,lower.tail,log.p) )
          }else{
            return( pchisq(-q/C,total_df,total_ncp,!lower.tail,log.p) )

          }
        }
      }

      # Although the cdf.func we return here does not make approximations in the tails, we
      # still need to extrapolation pieces below because when we integrate, we get much more stable
      # bounds by using our analytic integrals than trying to integrate over pchisq numerically.

      a.l <- a.r <- b.l <- b.r <- NULL

      if(f.eta[1]>0){
        # The following equations are found by simply matching the form of the extrapolated bound with the value and derivative
        # of the true cdf at the extrapolation point...this extrapolation helps make integrating the concentration inequality more stable
        e.r <- C*qchisq(-12*log(10),total_df,total_ncp,lower.tail = F,log.p=T)
        b.r <- exp(dchisq(e.r/C,total_df,total_ncp,log = T)-pchisq(e.r/C,total_df,total_ncp,F,T)-log(C))
        a.r <- -pchisq(e.r/C,total_df,total_ncp,F,T)-b.r*e.r
        e.l <- C*qchisq(-12*log(10),total_df,total_ncp,lower.tail = T,log.p=T)
        b.l <- -exp(dchisq(e.l/C,total_df,total_ncp,log = T)-pchisq(e.l/C,total_df,total_ncp,T,T)-log(C))
        a.l <- -pchisq(e.l/C,total_df,total_ncp,T,T)-b.l*e.l
      }

      if(f.eta[1]<0){
        # To get the below equations from the above, consider just negating x and taking 1- the form of the CDF to flip the CDF above twice.
        e.l <- - C*qchisq(-12*log(10),total_df,sum(delta2),lower.tail = F,log.p=T)
        b.l <- -exp(dchisq(-e.l/C,total_df,total_ncp,log = T)-pchisq(-e.l/C,total_df,total_ncp,F,T)-log(C))
        a.l <- -pchisq(-e.l/C,total_df,total_ncp,F,T)-b.l*e.l
        e.r <- -C*qchisq(-12*log(10),total_df,total_ncp,lower.tail = T,log.p=T)
        b.r <- -exp(dchisq(-e.r/C,total_df,total_ncp,log = T)-pchisq(-e.r/C,total_df,total_ncp,T,T)-log(C))
        a.r <- -pchisq(-e.r/C,total_df,total_ncp,T,T)+b.r*e.r
      }


      attr(cdf.func,"fft_used") <- FALSE
      attr(cdf.func,"f.eta") <- f.eta
      attr(cdf.func,"delta2") <- delta2
      attr(cdf.func,"df") <- df
      attr(cdf.func,"mu") <- sum(f.eta*(df+delta2))
      attr(cdf.func,"Q.sd") <- sum(2*(df+2*delta2)*f.eta^2)
      attr(cdf.func,"tail.features") <- list("support" = ifelse(f.eta[1] > 0,"pos.reals","neg.reals"),
                                             "extrapolation.point.l" = e.l,
                                             "extrapolation.point.r" = e.r,
                                             "a.l" = a.l,
                                             "b.l" = b.l,
                                             "a.r" = a.r,
                                             "b.r" = b.r)
      class(cdf.func) <- c("QFGaussCDF",class(cdf.func))
    return(cdf.func)
  }

  ncps <- delta2

  cdf <- calc.QFcdf(evals = f.eta, ncps = ncps, df = df, sigma = sigma, n = n, parallel.sapply = parallel.sapply)
  cdf.func <- wrap.QFcdf(cdf)

  attr(cdf.func,"fft_used") <- TRUE
  attr(cdf.func,"f.eta") <- f.eta
  attr(cdf.func,"delta2") <- delta2
  attr(cdf.func,"df") <- df
  attr(cdf.func,"sigma") <- sigma
  attr(cdf.func,"mu") <- cdf$mu
  attr(cdf.func,"Q.sd") <- cdf$Q.sd
  attr(cdf.func,"tail.features") <- list("support" = switch(cdf$type, mixed = "all.reals", pos = "pos.reals", neg = "neg.reals"),
                                         "extrapolation.point.l" = cdf$x[1],
                                         "extrapolation.point.r" = cdf$x[cdf$n],
                                         "a.l" = cdf$a.l,
                                         "b.l" = cdf$b.l,
                                         "a.r" = cdf$a.r,
                                         "b.r" = cdf$b.r)

  class(cdf.func) <- c("QFGaussCDF",class(cdf.func))
  cdf.func
}



#' Plot method for a QFGaussCDF object
#'
#' Plots the CDF computed by QFGauss.
#'
#' @param cdf a QFGaussCDF
#' @param ... additional parameters
#' @return There is nothing returned.
#'
#' @export
plot.QFGaussCDF <- function(x,...){
  if(class(x)[1]!="QFGaussCDF"){stop("x must be of class QFGaussCDF")}

  if(any(sapply(attr(x,"tail.features"),is.na))){warning("plotting domain truncated because at least one tail is missing: tail extrapolation in QFGauss failed, see ?QFGauss for details.")}

  xx <- calc.plotting.grid(x)

  old.par <- par(no.readonly = T)
  plot(xx, x(xx), type = "l", lwd=1.5, ylab = expression(CDF),xlab=expression(T[f])) # plot lower tail of CDF
  par(old.par)
}


#' Test function for a QFGaussCDF object
#'
#' Compares the CDF inferred by QFGauss to an empirical CDF.
#'
#' Four plots are produced.  The top-left plot overlays the CDF computed by QFGauss (in black) and the empirical CDF (in red) based on 10,000 samples.
#' The top-right plot shows the distance between the empirical and QFGauss-computed CDF and the corresponding ks.test p-value (two-sided alternative).
#' The ks.test p-value will be \code{NA} if \code{cdf} is missing one of its tails (see \code{\link{QFGauss}} for details).
#' The two bottom plots allow comparison of the empirical CDF (in red) with the computed CDF (in black) in each tail.
#'
#' @param cdf a QFGaussCDF
#' @param n.samps number of draws from the target distribution with which to construct the empirical CDF
#' @return Nothing is returned.
#' @seealso \code{\link{QFGauss}}, \code{\link{TestQFGaussBounds}}
#' @examples
#' TestQFGauss(QFGauss(c(1,3,4,-3)))
#'
#' @export
TestQFGauss <- function(cdf, n.samps = 1e4){

  if(class(cdf)[1]!="QFGaussCDF"){stop("cdf must be of class QFGaussCDF")}

  missing.tail <- FALSE

  if(any(is.na(c(attr(cdf,"tail.features")$a.r, attr(cdf,"tail.features")$b.r)))){
    missing.tail <- TRUE
    warning("testing domain truncated because cdf is missing its right tail: tail extrapolation in QFGauss failed, see ?QFGauss for details.")
  }

  if(any(is.na(c(attr(cdf,"tail.features")$a.l, attr(cdf,"tail.features")$b.l)))){
    missing.tail <- TRUE
    warning("testing domain truncated because cdf is missing its left tail: tail extrapolation in QFGauss failed, see ?QFGauss for details.")
  }

  f.eta <- attr(cdf,"f.eta")
  delta2 <- attr(cdf,"delta2")
  df <- attr(cdf,"df")
  sigma <- attr(cdf,"sigma")

  old.par <- par(no.readonly = T)
  par(mfrow=c(2,2))

  # break down samps to be sampled in chunks of 10,000
  if(n.samps > 1e4){
    n.samps.v <- diff(floor(seq(0,n.samps,length=ceiling(n.samps/1e4) )))
  } else {
    n.samps.v <- n.samps
  }
  samps <- as.list(1:length(n.samps.v))
  for(j in 1:length(n.samps.v)){
    samps[[j]] <- c(colSums(f.eta * matrix(rchisq(n.samps.v[j] * length(f.eta),df = df,ncp=delta2), nrow = length(f.eta)))) + rnorm(n.samps.v[j],sd = sigma)
  }
  samps <- do.call(c,samps)
  qecdf <- ecdf(samps)

  x <- calc.plotting.grid(cdf,range(samps))

  plot(x, cdf(x), type = "l", lwd=1.5, ylab = expression(CDF),xlab=expression(T[f]), main=expression(CDF~of~T[f]),las=1,bty="n") # plot lower tail of CDF
  lines(x, qecdf(x), lwd=1,col="red") # plot lower tail of CDF
  #legend("topleft",legend=c("QForm CDF","ECDF"),col=c("black","blue"),lty=c(1,3))

  if(missing.tail){ ks.pvalue <- NA }else{ ks.pvalue <- signif(ks.test(samps,cdf)$p.value,2) }

  plot(x,cdf(x)-qecdf(x), type="l",ylab = expression(CDF - ECDF), xlab=expression(T[f]), main=paste("Comparison to ECDF: ks.test p-value =", ks.pvalue),lwd=1.5,font.main=1,las=1,bty="n")
  abline(h=0)

  plot(x, -cdf(x, log.p = TRUE)/log(10), type = "l",
       ylab = expression(-log[10](CDF)), main=expression(-log[10](CDF)), xlab=expression(T[f]),las=1,bty="n")
  lines(x,-log10(qecdf(x)),col="red")
  #legend("topright",legend=c("QForm CDF","ECDF"),col=c("black","blue"),lty=c(1,3))


  plot(x, -cdf(x, lower.tail = FALSE, log.p = TRUE)/log(10), type = "l",
       ylab = expression(-log[10](1 - CDF)),main=expression(-log[10](1-CDF)), xlab=expression(T[f]),las=1,bty="n") # plot upper tail of CDF
  lines(x,-log1p(-qecdf(x))/log(10),col="red")
  #legend("topleft",legend=c("QForm CDF","ECDF"),col=c("black","blue"),lty=c(1,3))

  par(old.par)
}







