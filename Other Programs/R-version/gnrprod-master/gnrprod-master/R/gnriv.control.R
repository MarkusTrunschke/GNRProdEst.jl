#' Control parameters in \code{gnriv}
#' @description Allows the user to modify convergence parameters of Gauss Newton algorithm used in the \code{\link[gnrprod]{gnriv}} function
#' 
#' @param degree_w degree of Markov process for persistent productivity. Defaults to 3.
#' @param degree_tau degree of expansion for constant of integration. Defaults to the degree of first stage.
#' @param method the method of optimization passed to \code{\link[stats]{optim}}. Defaults to "BFGS." See \code{\link[stats]{optim}} under 'Details' for listing of available methods.
#' @param ... additional optional control parameters passed to \code{\link[stats]{optim}}. See \code{\link[stats]{optim}} for available parameters.
#' @return a list containing \code{degree} and \code{method} and any additional parameters in \code{...}. 
#' 
#' @usage gnriv.control(degree_w = 3, degree_tau = -1, method = "BFGS", ...)
#' 
#' @export

gnriv.control <- function(degree_w = 3, degree_tau = -1, method = "BFGS", ...) {
  return(list(degree_w = degree_w, degree_tau = degree_tau,
              method = method, ...))
}



