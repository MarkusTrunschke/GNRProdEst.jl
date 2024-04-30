#' Printing gross output function estimates
#' @description Print estimates of the parameters in a gross output function and names of the output, input, and data from \code{\link[gnrprod]{gnrprod}}.
#' 
#' @param x an object of class 'gnr'.
#' @param digits the number of significant figures to use for printing.
#' @param ... currently not used.
#' @return \code{print.gnr} has no return value and only prints a brief overview of elements contained in an object of class 'gnr' as described in the description.
#' @usage 
#' \method{print}{gnr}(x, digits = max(3L, getOption("digits") - 3L), ...)
#' @export
print.gnr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Gross Output Function:\n")
  cat("  output: ", deparse(x$call$output), "\n", sep = " ")
  cat("  fixed inputs: ", deparse(x$call$fixed), "\n", sep = " ")
  cat("  flexible inputs: ", deparse(x$call$flex), "\n", sep = " ")
  if (!is.null(x$call$data)) {
    cat("  data: ", deparse(x$call$data), "\n", sep = "")
  }
  cat("\n")
  cat("Estimates:\n")
  print(x$estimates$elas, digits = digits)
  invisible(x)
}


#' Summarizing gross output function estimates
#' @description Return a summary of the estimation routine for gross output functions from \code{\link[gnrprod]{gnrprod}}.
#' 
#' @param object an object of class 'gnr'.
#' @param ... currently not used.
#' @return a list of class 'summary.gnr' containing 14 elements:
#' \itemize{
#'  \item{\code{output_name}}{: the name of the output variable.}
#'  \item{\code{fixed_names}}{: a vector of the names of fixed input variables.}
#'  \item{\code{flex_name}}{: the name of the flexible input variable.}
#'  \item{\code{data}}{: `data` returned by \code{\link[gnrprod]{gnrprod}}.}
#'  \item{\code{data_name}}{: the name of `data`.}
#'  \item{\code{fs_conv}}{: a boolean indicating if convergence was achieved in the first stage.}
#'  \item{\code{ss_conv}}{: the convergence code of \code{\link[stats]{optim}} used in the second stage.}
#'  \item{\code{productivity}}{: matrix of the estimated total productivity.}
#'  \item{\code{fs_iter}}{: the number of iterations in the first stage.}
#'  \item{\code{fs_SSR}}{: sum of squared residuals in the first stage.}
#'  \item{\code{ss_iter}}{: the number of iterations in the second stage.}
#'  \item{\code{ss_val}}{: the value of the objective function in the second stage.}
#'  \item{\code{ss_iter}}{: the number of iterations in the second stage.}
#'  \item{\code{ss_mes}}{: the convergence message in the second stage.}
#' }
#' 
#' @usage 
#' \method{summary}{gnr}(object, ...)
#' @export
summary.gnr <- function(object, ...) {
  if (!inherits(object, "gnr")) {
    stop("'summary.gnr' called on an object not of class 'gnr'")
  }
  
  data <- object$data
  nfixed <- length(object$estimates$elas) - 1
  fixed_names <- colnames(data)[2:(nfixed + 1)]
  flex_name <- colnames(data)[nfixed + 2]
  output_name <- colnames(data)[1]
  
  if (is.null(object$call$data)) {
    data_name <- NULL 
  } else {
    data_name <- deparse(object$call$data)
  }
  
  fs_conv <- object$first_stage$convergence
  if (object$second_stage$optim_info$convergence == 0) {
    ss_conv <- TRUE
  } else {
    ss_conv <- FALSE
  }
  
  std_errors <- object$estimates$std_errors
  if (!is.null(std_errors)) {
    estimates <- cbind(object$estimates$elas, std_errors)
    colnames(estimates) <- c("Estimate", "Std. Error")
  } else {
    estimates <- as.matrix(object$estimates$elas)
    colnames(estimates) <- c("Estimate")
  }
  
  productivity <- object$data$productivity
  fs_iter <- object$first_stage$iterations
  fs_SSR <- object$first_stage$SSR
  ss_iter <- object$second_stage$optim_info$counts[1]
  ss_val <- object$second_stage$optim_info$value
  ss_mes <- object$second_stage$optim_info$message
  
  summary <- list(output_name = output_name,
                  fixed_names = fixed_names,
                  flex_name = flex_name,
                  data = data,
                  data_name = data_name,
                  fs_conv = fs_conv,
                  ss_conv = ss_conv,
                  estimates = estimates,
                  productivity = productivity,
                  fs_iter = fs_iter,
                  fs_SSR = fs_SSR,
                  ss_iter = ss_iter,
                  ss_val = ss_val,
                  ss_mes = ss_mes)
  class(summary) <- "summary.gnr"
  return(summary)
}


#' Printing a summary of gross output function estimation
#' @description Print a summary of the gross output function estimation routine from \code{\link[gnrprod]{gnrprod}}: names of output and inputs, summary statistics of the estimated productivity, the function estimates and standard errors if applicable, and convergence information.
#' 
#' @param x an object of class 'summary.gnr'.
#' @param digits the number of significant figures to use for printing.
#' @param ... currently not used.
#' 
#' @return \code{print.gnr} has no return value and only prints the elements contained in an object of class 'summary.gnr' as described in the description.
#' 
#' @usage 
#' \method{print}{summary.gnr}(x, digits = max(3L, getOption("digits") - 3L), ...)
#' @export
print.summary.gnr <- function(x, digits = max(3L, getOption("digits") - 3L),
                              ...) {
  cat("Gross Output Function:\n")
  cat("  output: ", deparse(x$output_name), "\n", sep = "")
  cat("  fixed inputs: ", paste(deparse(x$fixed_names), sep = "", collapse = ", "), "\n")
  cat("  flexible inputs: ", deparse(x$flex_name), "\n", sep = "")
  if (!is.null(x$data_name)) {
    cat("  data: ", x$data_name, "\n", sep = "")
  }
  cat("\n")
  cat("Total Productivity:\n")
  print(summary(x$productivity))
  cat("\n")
  cat("Estimates:\n")
  stats::printCoefmat(x$estimates, digits = digits)
  cat("\n")
  cat("First-Stage Convergence:", x$fs_conv, "with", x$fs_iter, "iterations. Sum of Squared Residuals:", x$fs_SSR, "\n", sep = " ")
  if (x$ss_conv) {
    cat("Second-Stage Convergence:", x$ss_conv, "with", x$ss_iter, "iterations. Value:", x$ss_val, "\n", sep = " ")
  } else {
    cat("Second-Stage Convergence:", FALSE, x$ss_mes, "\n", sep = " ")
  }
}

#' Print gross output function estimates
#' @description Print or return a numeric matrix of the estimated parameters from an object of class 'gnr'.
#' 
#' @param object an object of class 'gnr'.
#' @param ... currently not used.
#' 
#' @return the named vector of parameter estimates contained in an object of class 'gnr'.
#' 
#' @usage 
#' \method{coef}{gnr}(object, ...)
#' @export
coef.gnr <- function(object, ...) object$estimates$elas


#' Printing first stage estimate
#' @description Print estimate of the flexible input elasticity, the sum of squared residuals, and convergence status from \code{\link[gnrprod]{gnrflex}}.
#' 
#' @param x an object of class 'gnrflex'.
#' @param digits the number of significant figures to use for printing.
#' @param ... currently not used.
#' 
#' @return \code{print.gnrflex} has no return value and only prints a brief overview of elements contained in an object of class 'gnrflex' as described in the description.
#' 
#' @usage
#' \method{print}{gnrflex}(x, digits = max(3L, getOption("digits") - 3L), ...)
#' @export
print.gnrflex <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Flexible Input Elasticity:\n")
  flex_elas_avg <- c(mean(x$elas$flex_elas))
  names(flex_elas_avg) <- colnames(x$elas$flex_elas)
  print(flex_elas_avg, digits = digits)
  cat("\n")
  cat("First-Stage Sum of Squared Residuals: ", x$elas$SSR, "\n", sep = "")
  cat("Convergence: ", x$elas$convergence, "\n", sep = "")
  invisible(x)
}

#' Printing second stage estimates
#' @description Print estimates of the fixed input elasticities, productivity, objective function value, and convergence status from \code{\link[gnrprod]{gnriv}}.
#' 
#' @param x an object of class 'gnriv'.
#' @param digits the number of significant figures to use for printing.
#' @param ... currently not used.
#' 
#' @return \code{print.gnriv} has no return value and only prints a brief overview of elements contained in an object of class 'gnriv' as described in the description.
#' 
#' @usage 
#' \method{print}{gnriv}(x, digits = max(3L, getOption("digits") - 3L), ...)
#' @export
print.gnriv <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Fixed Input Elasticity:\n")
  fixed_elas_avg <- apply(x$fixed_elas, 2, FUN = mean)
  print(fixed_elas_avg, digits = digits)
  cat("\n")
  cat("Total Productivity:\n")
  print(summary(x$productivity))
  cat("\n")
  cat("Second-Stage Objective Function Value: ", x$optim_info$value, "\n", sep = "")
  if (x$optim_info$convergence == 0) {
    cat("Convergence: ", TRUE, "\n", sep = " ")
  } else {
    cat("Convergence:", FALSE, x$optim_info$message, "\n", sep = " ")
  }
  invisible(x)
}








