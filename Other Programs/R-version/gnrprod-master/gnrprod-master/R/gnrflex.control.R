#' Control parameters in \code{gnrflex}
#' @description Allows the user to modify convergence parameters of Gauss Newton algorithm used in the \code{\link[gnrprod]{gnrflex}} function.
#'
#' @param degree degree of share regression polynomial. Defaults to 3.
#' @param maxit maximum number of iterations. Defaults to 100.
#' @param reltol relative convergence tolerance. Defaults to 1e-5.
#' @param initial_step a scaling parameter specifying the initial step-size factor used in each iteration of the Gauss-Newton algorithm. \code{initial_step} is halved in each convergence step.
#' @param min_factor the minimum value that the step-size factor can take on in the convergence step of any iteration of the Gauss-Newton algorithm.
#' @return a list containing five elements: \code{degree}, \code{maxit}, \code{reltol}, \code{initial_step}, and \code{min_factor}.
#' 
#' @usage gnrflex.control(degree = 3, maxit = 100, reltol = 1e-5,
#'                        initial_step = 100, min_factor = 1e-5)
#' 
#' @export

gnrflex.control <- function(degree = 3, maxit = 100, reltol = 1e-5,
                            initial_step = 100, min_factor = 1e-5) {
  if (maxit <= 0 || !is.numeric(maxit)) {
    stop("maximum iteration count must be integer > 0")
  }

  if (reltol <= 0 || !is.numeric(reltol)) {
    stop("relative tolerance must be numeric > 0")
  }

  if (initial_step <= 0 || !is.numeric(initial_step)) {
    stop("initial_step must be numeric > 0")
  }

  if (min_factor < 0 || !is.numeric(min_factor)) {
    stop("min_factor must be numeric >= 0")
  }

  list(degree = degree, maxit = maxit, reltol = reltol,
       initial_step = initial_step, min_factor = min_factor)
}
