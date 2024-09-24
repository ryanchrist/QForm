#' QForm: A package for fast, safe CDF/PDF estimation for generalized chi-square random varaibles and screening with p-value bounds for quadratic forms.
#'
#' QForm returns the CDF and PDF for the generalized chi-squared distribution with numerical accuracy deep into the tail (see \code{\link{QFGauss}}). It can also provide reliable upper and lower bounds on the CDF when only some of the chi-square terms are known (see \code{\link{QFGaussBounds}}).
#' By using the fast Fourier transform in combination with novel concentration inequalities and various adjustments for numerical precision, QForm function is faster and more reliable than Davie's method and related approaches, especially when the returned CDF or PDF is to be evaluated at many points.
#'
#' Initially motivated by genome-wide association studies (GWAS), QForm is aimed at obtaining upper and lower bounds on p-values for test statistics with a limiting distribution of the form
#' of \eqn{Q_f = T_f + R_f} where only the CDF of \eqn{T_f} is known. These random variables have the form
#'
#' \deqn{T_f = \sum\limits_{i \in \mathcal{T}} f\left(\eta_i \right) A_i + \sigma Z_0}{T_f = \Sigma_{i \in T} f (\eta_i) A_i + \sigma Z_0},
#' \deqn{R_f = \sum\limits_{i \in \mathcal{R}} f\left(\eta_i \right) A_i}{T_f = \Sigma_{i \in R} f (\eta_i) A_i},
#' where each \eqn{A_i \sim \chi^2_{a_i}\left(\delta^2_i\right)}{A_i \sim \chi^2_{a_i}\left(\delta^2_i\right)} and \eqn{Z_0 \sim N(0,1)}{Z_0 \sim N(0,1)}, all mututally independent,
#' and \eqn{a_i = 1}{a_i = 1} for all \eqn{i \in \mathcal{R}}{i \in \mathcal{R}}. We aim to remove this final restriction in future work.
#'
#' In the genomics literature, SKAT and related methods have limiting distributions of this form.  In the machine learning and kernel methods literature, other popular test statistics share this limiting distribution, among them the Hilbert-Schimidt Information Criterion (HSIC).
#' Approximate methods have emerged in the genomics (eg: FastSKAT) and kernel methods literature, have emerged based on the idea of using a top-k SVD to obtain \eqn{T_f} and then attempt to approximate the contribution from \eqn{R_f} using a single random variable that matches some of the moments of \eqn{R_f}.
#' However, we've found that in several applications, such approximations of \eqn{R_f} can lead to p-value estimates that are off by orders of magnitude.  We take a concentration-inequality based approach to bounding the potential contribution of \eqn{R_f}
#' to the overall distribution of \eqn{Q_f}, allowing us to obtain exact upper and lower bounds on the p-value that can allow users to quickly discard observations (eg: genomic loci) that could never be significant while concentrating further computational resources
#' on more precisely evaluating the p-value at loci that could still potentially be interesting/significant.
#'
#' Our implementation features two main new functions.  First, we do not rely on \code{CompQuadForm}, which implements Davie's method but as such has difficult-to-tune parameters and can often fail for pvalues smaller than 1e-16.
#' Davie's method is based on a more general integral transform that relates the CDF of a random variable to its characteristic function, but predates the fast fourier transform.  We make use of the same identity as Davies, but by combining it with the FFT, obtain the CDF of random variables of the form of \eqn{Q_f}
#' at many points in parallel (implemented in \code{\link{QFGauss}}).
#'
#' Given the CDF produced by \code{\link{QFGauss}}, we apply a set of analytic and numerical intergration routines to \eqn{T_f} to calculate our p-value bounds for \eqn{Q_f} (implemented in \code{\link{QFGaussBounds}}).
#'
#'
#' @section QForm functions:
#' \code{\link{QFGauss}}
#' \code{\link{TestQFGauss}}
#' \code{\link{QFGaussBounds}}
#' \code{\link{TestQFGaussBounds}}
#'
#' @aliases QForm-package
#' @name QForm
NULL


