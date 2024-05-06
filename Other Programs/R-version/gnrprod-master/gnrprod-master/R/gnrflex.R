#' Estimate flexible input elasticity: Gandhi, Navarro, Rivers (GNR) share regression; first stage
#' @description The \code{gnrflex} function implements the first stage (share
#' regression) of the GNR production function estimation routine,
#' nonparametrically identifying the flexible input elasticity of the
#' production function. This function is called within the main wrapper
#' function \code{\link[gnrprod]{gnrprod}}. If the production-related inputs
#' are characters, a \code{\link[base]{data.frame}} must be specified under
#' \code{data}. Alternatively, matrices/vectors may be directly specified
#' without specifying \code{data}. \code{gnrprod} currently supports only one
#' flexible input. The parameters are optimized using the Gauss-Newton
#' algorithm. \code{gnrflex} currently supports only one flexible input.
#'
#' For details, see Gandhi, Navarro, and Rivers (2020).
#'
#' @param output name (character) of variable of log gross output in data or a numeric vector.
#' @param fixed name (character or character vector) of variables of log fixed inputs in data or a numeric matrix.
#' @param flex name (character) of variable of log flexible input in data or a numeric vector.
#' @param share name (character) of variable of log intermediate input's revenue share in data or a numeric vector.
#' @param id name (character) of variable of firm id in data or a numeric vector.
#' @param time name (character) of variable of time in data or a numeric vector.
#' @param data \code{\link[base]{data.frame}} containing all variables with names specified by arguments above (left empty if arguments above are vector/matrix rather than strings).
#' @param control an optional list of convergence settings. See \code{\link[gnrprod]{gnrflex.control}} for listing.
#' @return a list of class 'gnrflex' containing three elements:
#'
#' \code{elas}: a list containing six elements describing the share regression:
#' \itemize{
#'  \item{\code{flex_elas}}{: a numeric vector of the estimated flexible input elasticity for each observation.}
#'  \item{\code{coefficients}}{: a numeric vector of the coefficients of the estimator scaled by a constant. See Gandhi, Navarro, and Rivers (2020, p. 2994, equation (21)).}
#'  \item{\code{residuals}}{: a numeric vector of the residuals.}
#'  \item{\code{SSR}}{: sum of squared residuals.}
#'  \item{\code{iterations}}{: number of iterations performed.}
#'  \item{\code{convergence}}{: boolean indicating whether convergence was achieved.}
#' }
#'
#' \code{arg}: a list containing eight elements to be passed to the second stage function \code{\link[gnrprod]{gnriv}}:
#' \itemize{
#'  \item{\code{input}}{: a numeric matrix (S3: \code{\link[stats]{poly}}) of the polynomial expansion of all inputs.}
#'  \item{\code{input_degree}}{: a numeric matrix corresponding to \code{input} denoting each vector's degree.}
#'  \item{\code{all_input}}{: a numeric matrix of the inputs without polynomial expansion.}
#'  \item{\code{big_Y}}{: a numeric vector of persistent productivity minus the constant of integration. See Gandhi, Navarro, and Rivers (2020, p. 2991, equation (16)).}
#'  \item{\code{D_coef}}{: a numeric vector equaling \code{coef} divided by an estimate of the constant.}
#'  \item{\code{id}}{: a numeric vector of the firm ids.}
#'  \item{\code{time}}{: a numeric vector of time.}
#'  \item{\code{degree}}{: the degree of the share regression.}
#'  \item{\code{fixed_names}}{: the names of fixed inputs. To be used in the second stage.}
#' }
#'
#' \code{control}: the list of convergence control parameters. See \code{\link[gnrprod]{gnrflex.control}} for available parameters.
#' 
#' @usage gnrflex(output, fixed, flex, share, id, time, data, control)
#' @examples 
#' require(gnrprod)
#' data <- colombian
#' industry_311_flex <- gnrflex(output = "RGO", fixed = c("L", "K"),
#'                              flex = "RI", share = "share", id = "id",
#'                              time = "year", data = data,
#'                              control = list(degree = 2, maxit = 200))
#' @references Gandhi, Amit, Salvador Navarro, and David Rivers. 2020. "On the Identification of Gross Output Production Functions." *Journal of Political Economy*, 128(8): 2973-3016. \doi{10.1086/707736}.
#' 
#' Davidson, Russell, James G. MacKinnon. 1993. "The Gauss-Newton Regression." In *Estimation and Inference in Econometrics*, 176-207. New York: Oxford University Press.
#' @export

gnrflex <- function(output, fixed, flex, share, id, time, data, control) {
  
  ctrl <- gnrflex.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  
  output <- get_matrix(output, data)
  fixed <- get_matrix(fixed, data)
  flex <- get_matrix(flex, data)
  share <- get_matrix(share, data)
  id <- get_matrix(id, data)
  time <- get_matrix(time, data)
  
  complete_obs <- stats::complete.cases(cbind(output, fixed, flex, id, time,
                                              share))
  if (sum(complete_obs) != length(output)) {
    output <- output[complete_obs, , drop = FALSE]
    fixed <- fixed[complete_obs, , drop = FALSE]
    flex <- flex[complete_obs, , drop = FALSE]
    id <- id[complete_obs, , drop = FALSE]
    time <- time[complete_obs, , drop = FALSE]
    share <- share[complete_obs, , drop = FALSE]
    warning("'output', 'fixed', 'flex', 'share', 'id', and 'time' contains missing values: observations omitted")
  }
  
  all_input <- cbind(fixed, flex)
  poly_input <- stats::poly(all_input, degree = ctrl$degree, raw = TRUE)
  input_degrees <- sapply(colnames(poly_input), FUN = function(x) {
    a <- base::strsplit(x, split = "[.]")[[1]]
  })
  input_degrees <- apply(input_degrees, 2, as.numeric)

  gamma_denom <- rbind(1, as.matrix(input_degrees[nrow(input_degrees), ] + 1))
  start_reg <- stats::lm(share ~ poly_input)

  constant <- start_reg$fitted.values - stats::coef(start_reg)[1]
  constant <- -min(constant, na.rm = TRUE) + 0.1

  start <- c(constant, (coef(start_reg)[-1]))
  share_reg <- gauss_newton_reg(start = start, data = poly_input,
                                share = share, control = ctrl)

  coef <- share_reg[[1]]
  i_elas <- log(pred_fs(coef, poly_input))
  errors <- i_elas - share
  colnames(errors) <- "fs_residuals"
  mean_exp_err <- mean(exp(errors))
  i_elas <- exp(i_elas - log(mean_exp_err))
  colnames(i_elas) <- colnames(flex)
  gamma <- as.matrix(coef / mean_exp_err)
  flex_gamma <- gamma / gamma_denom

  integ_G_I <- pred_fs(flex_gamma, poly_input)
  integ_G_I <- integ_G_I * flex
  big_Y <- as.matrix(output - errors - integ_G_I)
  colnames(big_Y) <- "big_Y"

  fs_elas <- list("flex_elas" = i_elas,
                  "coefficients" = coef,
                  "residuals" = errors,
                  "SSR" = c(share_reg$SSR),
                  "iterations" = share_reg$iterations,
                  "convergence" = share_reg$convergence)

  fs_arg <- list("input" = poly_input,
                 "input_degree" = input_degrees,
                 "all_input" = all_input,
                 "big_Y" = big_Y,
                 "D_coef" = gamma[-1, ] / gamma_denom[-1, ],
                 "id" = id,
                 "time" = time,
                 "degree" = ctrl$degree,
                 "fixed_names" = colnames(fixed))

  fs_return <- list("elas" = fs_elas, "arg" = fs_arg, "control" = ctrl)
  class(fs_return) <- "gnrflex"
  return(fs_return)
}

gauss_newton_reg <- function(start, data, share, control) {
  iter <- 0
  call_start <- start

  inputs_1 <- data.frame(rep(1, nrow(data)), data)
  names(inputs_1)[1] <- "constant"

  while (iter < control$maxit) {
    initial_pred <- pred_fs(call_start, data)

    X <- as.matrix(inputs_1 / initial_pred)

    initial_errors <- as.matrix(cbind(share - log(initial_pred)))
    initial_SSR <- t(initial_errors) %*% initial_errors

    new_start <- call_start + (solve(t(X) %*% X) %*% t(X) %*% initial_errors)
    new_pred <- pred_fs(new_start, data)
    suppressWarnings(new_errors <- cbind(share - log(new_pred)))
    new_SSR <- t(new_errors) %*% new_errors

    initial_step <- control$initial_step
    min_factor <- control$min_factor

    while ((is.na(new_SSR) || new_SSR > initial_SSR)
           & initial_step >= min_factor) {
      initial_step <- initial_step / 2
      new_start <- call_start +
        initial_step * (solve(t(X) %*% X) %*% t(X) %*% initial_errors)
      new_pred <- pred_fs(new_start, data)
      suppressWarnings(new_errors <- cbind(share - log(new_pred)))
      new_SSR <- t(new_errors) %*% new_errors
    }

    conv_bool <- TRUE
    for (i in 1:length(new_start)) {
      if (initial_step * abs(new_start[i]) > control$reltol *
          (abs(new_start[i]) + 1e-3)) {
        conv_bool <- FALSE
        break
      }
    }

    iter <- iter + 1
    if (conv_bool) {
      return_list <- list(new_start, iter, new_SSR, conv_bool, control)
      names(return_list) <- c("share_reg_coef", "iterations", "SSR",
                              "convergence")
      return(return_list)
    }
    call_start <- new_start
  }
  warning("share regression failed to converge")
  return_list <- list(new_start, iter, new_SSR, conv_bool)
  names(return_list) <- c("share_reg_coef", "iterations", "SSR", "convergence")
  return(return_list)
}


pred_fs <- function(start, data) {
  matrix <- as.matrix(cbind(rep(1, nrow(data)), data))
  new_m <- matrix %*% start
  return(new_m)
}



