#' Estimate production functions and productivity: Gandhi, Navarro, and Rivers (2020)
#' @description The \code{gnrprod} function is the front end of the
#' \code{gnrprod} package. It estimates production functions and productivity
#' in two stages: \code{\link[gnrprod]{gnrflex}} (estimate flexible input elasticity) and
#' \code{\link[gnrprod]{gnriv}} (estimate fixed input elasticities and productivity).
#' If the production-related inputs are characters, a \code{\link[base]{data.frame}}
#' must be specified under \code{data}. Alternatively, matrices/vectors may be
#' directly specified without specifying \code{data}. \code{gnrprod} currently
#' supports only one flexible input.
#'
#' @param output name (character) of variable of log gross output in data or a numeric vector.
#' @param fixed name (character or character vector) of variables of log fixed inputs in data or a numeric matrix.
#' @param flex name (character) of variable of log flexible input in data or a numeric vector.
#' @param share name (character) of variable of log intermediate input's revenue share in data or a numeric vector.
#' @param in_price optional (required if \code{share} is not specified) name (character) of variable of common flexible input price or a numeric vector.
#' @param out_price optional (required if \code{share} is not specified) name (character) of variable of common output price or a numeric vector.
#' @param id name (character) of variable of firm ID in data or a numeric vector.
#' @param time name (character) of variable of time in data or a numeric vector.
#' @param data \code{\link[base]{data.frame}} containing all variables with names specified by arguments above (left empty if arguments above are vector/matrix rather than strings).
#' @param B number of bootstrap repetitions to retrieve standard errors of elasticity estimates. By default, \code{gnrprod} does not bootstrap, i.e., \code{B = NULL}. Setting \code{B > 1} will output bootstrapped standard errors. 
#' @param fs_control an optional list of convergence settings of the first stage. See \code{\link[gnrprod]{gnrflex.control}} for listing.
#' @param ss_control an optional list of convergence settings of the second stage. See \code{\link[gnrprod]{gnriv.control}} for listing.
#' @param ... additional optional arguments to be passed to \code{\link[stats]{optim}} in the second stage.
#' @return a list of class 'gnr' with five elements:
#' 
#' \code{estimates}: a list with two elements: \code{elas} the parameter estimates and \code{std_errors} the standard errors.
#'
#' \code{data}: a \code{\link[base]{data.frame}} containing: \code{output}, \code{fixed}, \code{flex}, \code{share}, \code{id}, \code{time}, estimated elasticities for each observation, estimated productivity, and first stage residuals.
#'
#' \code{first_stage}: a list containing five elements describing the share regression (first stage):
#' \itemize{
#'  \item{\code{coefficients}}{: a numeric vector of the coefficients of the first stage estimator scaled by a constant. See Gandhi, Navarro, and Rivers (2020, p. 1994, equation (21)).}
#'  \item{\code{SSR}}{: sum of squared residual.}
#'  \item{\code{iterations}}{: number of iterations performed.}
#'  \item{\code{convergence}}{: boolean indicating whether convergence was achieved.}
#'  \item{\code{control}}{: list of convergence control parameters (see \code{\link[gnrprod]{gnrflex.control}}).}
#' }
#'
#' \code{second_stage}: a list containing four elements describing the second stage:
#' \itemize{
#'  \item{\code{optim_method}}{: the method for optimization. Defaults to 'BFGS'. See \code{\link[stats]{optim}} for a listing of available methods.}
#'  \item{\code{optim_info}}{: the returned list of the \code{\link[stats]{optim}} function estimating the coefficients of the constant of integration. See Gandhi, Navarro, and Rivers (2020, p. 1994, equation (21)).}
#'  \item{\code{optim_control}}{: the list of control parameters passed to \code{\link[stats]{optim}}.}
#'  \item{\code{degree_w}}{: degree of Markov process for persistent productivity.}
#'  \item{\code{degree_tau}}{: degree of expansion for constant of integration.}
#' }
#' 
#' \code{call}: the function call.
#'
#' @usage gnrprod(output, fixed, flex, share, in_price = NULL,
#'                out_price = NULL, id, time, data, B = NULL,
#'                fs_control = NULL, ss_control = NULL, ...)
#' 
#' @examples
#' require(gnrprod)
#' data <- colombian
#' industry_311 <- gnrprod(output = "RGO", fixed = c("L", "K"),
#'                         flex = "RI", share = "share", id = "id",
#'                         time = "year", data = data,
#'                         fs_control = list(degree = 2, maxit = 200),
#'                         ss_control = list(trace = 1))
#'                         
#' @references Gandhi, Amit, Salvador Navarro, and David Rivers. 2020. "On the Identification of Gross Output Production Functions." *Journal of Political Economy*, 128(8): 2973-3016. \doi{10.1086/707736}.
#' @export

gnrprod <- function(output, fixed, flex, share, in_price = NULL,
                    out_price = NULL, id, time, data, B = NULL,
                    fs_control = NULL, ss_control = NULL, ...) {

  cl <- match.call()
  
  output <- get_matrix(output, data)
  fixed <- get_matrix(fixed, data)
  flex <- get_matrix(flex, data)
  id <- get_matrix(id, data)
  time <- get_matrix(time, data)
  
  if (!missing(share)) {
    share <- get_matrix(share, data)
  } else if (!missing(in_price) && !missing(out_price)) {
    share <- (in_price + flex) - (out_price + output)
    in_price <- get_matrix(in_price, data)
    out_price <- get_matrix(out_price, data)
    colnames(share) <- "share1"
  } else {
    stop("must specify either share or both intermediate-input price and output price")
  }

  complete_obs <- stats::complete.cases(cbind(output, fixed, flex, share, id,
                                              time))
  
  if (sum(complete_obs) != length(output)) {
    output <- output[complete_obs, , drop = FALSE]
    fixed <- fixed[complete_obs, , drop = FALSE]
    flex <- flex[complete_obs, , drop = FALSE]
    id <- id[complete_obs, , drop = FALSE]
    time <- time[complete_obs, , drop = FALSE]
    share <- share[complete_obs, , drop = FALSE]
    warning("'output', 'fixed', 'flex', 'share', 'id', and 'time' contains incomplete observations: observations omitted")
  }
  
  if (!is.null(in_price) && !is.null(out_price)) {
    in_price <- in_price[complete_obs, , drop = FALSE]
    out_price <- out_price[complete_obs, , drop = FALSE]
    mf <- cbind(output, fixed, flex, in_price, out_price, share,
                           id, time)
  } else {
    mf <- cbind(output, fixed, flex, share, id, time)
  }
  
  mf <- mf[order(mf[, ncol(mf) - 1], mf[, ncol(mf)]), , drop = FALSE]
  
  
  gnr_flex <- gnrflex(output = output, fixed = fixed, flex = flex,
                      share = share, id = id, time = time,
                      control = fs_control)

  gnr_iv <- gnriv(object = gnr_flex, control = ss_control, ...)

  boot_sd <- NULL
  if (!missing(B) && B > 1) {
    id_unique <- unique(id)
    boot_elas <- lapply(1:B, FUN = function(i) {
      boot_ids <- sample(id_unique, length(id_unique), replace = TRUE)
      boot_df <- do.call(rbind, lapply(boot_ids, function (x) {
        mf[mf[, colnames(id)] == x,]
      }))
      
      boot_output <- boot_df[, 1]
      boot_fixed <- boot_df[, 2:(ncol(fixed) + 1)]
      boot_flex <- boot_df[, (ncol(fixed) + 2):(ncol(fixed) + 1 + ncol(flex))]
      boot_share <- boot_df[, ncol(boot_df) - 2]
      boot_id <- boot_df[, ncol(boot_df) - 1]
      boot_time <- boot_df[, ncol(boot_df)]
      
      boot_fs <- suppressWarnings(gnrflex(output = boot_output,
                                         fixed = boot_fixed, flex = boot_flex,
                                         share = boot_share, id = boot_id,
                                         time = boot_time,
                                         control = fs_control))
      
      flex_elas <- mean(boot_fs$elas$flex_elas)
      
      boot_ss <- suppressWarnings(gnriv(object = boot_fs,
                                        control = ss_control, ...))
      
      fixed_elas <- apply(boot_ss$fixed_elas, 2, mean)
      
      return(c(fixed_elas, flex_elas))
    })
    boot_est <- do.call(cbind, boot_elas)
    boot_sd <- apply(boot_est, 1, stats::sd)
  }
  
  fixed_elas <- gnr_iv$fixed_elas
  flex_elas <- gnr_flex$elas$flex_elas
  elas <- data.frame(cbind(fixed_elas, flex_elas))
  input_names <- c(colnames(fixed), colnames(flex))
  colnames(elas) <- input_names
  mf <- cbind(mf, elas, gnr_iv$productivity, gnr_flex$elas$residuals)
  colnames(mf)[(ncol(mf) - 1):ncol(mf)] <- c("productivity", "flex_resid")

  fs_return <- list("coefficients" = gnr_flex$elas$coefficients,
                    "SSR" = gnr_flex$elas$SSR,
                    "iterations" = gnr_flex$elas$iterations,
                    "convergence" = gnr_flex$elas$convergence,
                    "control" = gnr_flex$control)

  ss_return <- list("optim_method" = gnr_iv$control$method,
                    "optim_info" = gnr_iv$optim_info,
                    "optim_control" = gnr_iv$control$optim_control,
                    "degree_constant" = gnr_iv$control$degree_tau,
                    "degree_productivity" = gnr_iv$control$degree_w)

  return_average_elas <- apply(elas, MARGIN = 2, FUN = mean)
  
  return_sd <- boot_sd
  if (!is.null(boot_sd)) {
    names(return_sd) <- paste(input_names, "se", sep = "_")
  }
  param = list("elas" = return_average_elas, "std_errors" = return_sd)
  
  gnr_out <- list(param, mf, fs_return, ss_return, cl)
  names(gnr_out) <- c("estimates", "data", "first_stage", "second_stage",
                      "call")
  class(gnr_out) <- "gnr"
  return(gnr_out)
}

get_matrix <- function(x, data) {
  if (is.null(x)) {
    return(NULL)
  }
  
  if (is.character(x) && missing(data)) {
    stop(paste(deparse(substitute(x)), "is a string but argument `data` is missing"))
  }
  
  if (is.character(x)) {
    col <- as.matrix(data[, x])
    colnames(col) <- x
  } else {
    col <- as.matrix(x)
    colnames(col) <- colnames(x, do.NULL = FALSE,
                              prefix = paste(substitute(x)))
  }
  return(col)
}






