#' Estimate fixed input elasticity and total productivity: Gandhi, Navarro, Rivers (GNR) lags as instruments; second stage
#' @description The \code{gnriv} function implements the second stage of the
#' GNR production function estimation routine, nonparametrically identifying
#' the fixed input elasticities of the production function and total
#' productivity. This function accepts an object of class 'gnrflex'. The
#' parameters are optimized using the function \code{\link[stats]{optim}}.
#'
#' For details, see Gandhi, Navarro, and Rivers (2020).
#'
#' @param object object of class 'gnrflex'.
#' @param control an optional list of convergence settings. See \code{\link[gnrprod]{gnriv.control}} for listing.
#' @param ... additional optional arguments passed to optim.
#' @return a list of class 'gnriv' containing three elements:
#'
#' \code{fixed_elas}: a numeric matrix of estimated elasticities of fixed inputs for each observation.
#'
#' \code{productivity}: a numeric vector of estimated total productivity.
#'
#' \code{control}: the list of convergence control parameters. See \code{\link[gnrprod]{gnriv.control}} for available parameters.
#' 
#' @usage gnriv(object, control, ...)
#' 
#' @examples 
#' require(gnrprod)
#' data <- colombian
#' \dontrun{
#' industry_311_flex <- gnrflex(output = "RGO", fixed = c("L", "K"),
#'                              flex = "RI", share = "share", id = "id",
#'                              time = "year", data = data,
#'                              control = list(degree_w = 2, maxit = 200))
#' 
#' industry_311_fixed <- gnriv(industry_311_flex,
#'                             control = list(trace = 1))}
#' @importFrom data.table "data.table"
#' @importFrom data.table ".SD"
#' @importFrom data.table "shift"
#' @references Gandhi, Amit, Salvador Navarro, and David Rivers. 2020. "On the Identification of Gross Output Production Functions." *Journal of Political Economy*, 128(8): 2973-3016. \doi{10.1086/707736}.
#' @export


gnriv <- function(object, control, ...) {
  # Check if object has correct class
  if (attr(object, "class") != "gnrflex") {
    stop("object must be of class gnrflex")
  }
  
  # Hand over controls if any
  ctrl <- gnriv.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
    if (length(ctrl) == 3) {
      optim.control <- NULL
    } else {
      optim.control <- ctrl[4:length(ctrl)]
    }
  } else {
    optim.control <- NULL
  }
  
  degree_w <- ctrl[[1]]
  degree_tau <- ctrl[[2]]
  method <- ctrl[[3]]
  
  if (degree_tau <= 0) {
    degree_tau <- object$control$degree
  }
  
  if (degree_tau == object$control$degree) {
    all_input <- object$arg$input
    orig_input <- all_input
    input_degree <- object$arg$input_degree
    orig_input_degree <- input_degree
  } else {
    all_input <- object$arg$all_input
    all_input <- stats::poly(all_input, degree = degree_tau, raw = TRUE)
    orig_input <- object$arg$input
    input_degree <- sapply(colnames(all_input), FUN = function(x) {
      a <- base::strsplit(x, split = "[.]")[[1]]
    })
    input_degree <- apply(input_degree, 2, as.numeric)
    orig_input_degree <- object$arg$input_degree
  }
  
  pred <- sapply(1:ncol(input_degree), FUN = function(i) {
    if (input_degree[nrow(input_degree), i] == 0) {
      return(all_input[, i, drop = FALSE])
    }
  })
  pred <- do.call(cbind, pred) # Get columns corresponding to polynomias from fixed inputs

  id <- object$arg$id
  time <- object$arg$time
  big_Y <- object$arg$big_Y
  errors <- object$elas$residuals

  # Get starting values
  constant_reg <- stats::lm(big_Y ~ as.matrix(pred)) # Regress weird_Y on fixed input polynomials. Not sure if he adds a constant. GNR do and I follow them for now.
  

  coefficients <- stats::coef(constant_reg)[2:(ncol(pred) + 1)]
  names(coefficients) <- colnames(pred)
  
  fixed_base <- data.table::data.table(id, time, pred, big_Y)
  colnames(fixed_base)[1:2] <- c("id", "time")
  fixed_base <- fixed_base[order(fixed_base$id, fixed_base$time), ]
  fixed_lag <- fixed_base[, data.table::shift(.SD, n = 1), by = id,
                          .SDcols = 3:ncol(fixed_base), drop = FALSE] # Get lags of the fixed inputs.
  
  complete_obs <- stats::complete.cases(cbind(fixed_base, fixed_lag)) # Drop those who have missings b/c of lags
  
  fixed_base <- as.matrix(fixed_base)[complete_obs, -1:-2, drop = FALSE]
  fixed_lag <- as.matrix(fixed_lag)[complete_obs, -1, drop = FALSE]
  
  big_Y_base <- fixed_base[, ncol(fixed_base)]
  fixed_base <- fixed_base[, -c(ncol(fixed_base))]
  big_Y_lag <- fixed_lag[, ncol(fixed_lag)]
  fixed_lag <- fixed_lag[, -c(ncol(fixed_lag))]

  # GMM starts here
  constant_gmm <- stats::optim(par = coefficients, fn = constant_moments,
                               data = fixed_base, big_Y_base = big_Y_base,
                               big_Y_lag = big_Y_lag, lag_data = fixed_lag,
                               degree = degree_w, method = method,
                               control = optim.control, ...)
  
  # Not sure why he does this. Probably to fill up the return vector if optim.control did not specify anything that deviates from default.
  opt_ctrl <- list(trace = 0, # Don't display progress
                   fnscale = 1, # An overall scaling to be applied to the value of fn and gr during optimization. If negative, turns the problem into a maximization problem. Optimization is performed on fn(par)/fnscale.
                   parscale = rep.int(1, length(coefficients)), # A vector of scaling values for the parameters. Optimization is performed on par/parscale and these should be comparable in the sense that a unit change in any element produces about a unit change in the scaled value. Not used (nor needed) for method = "Brent".
                   ndeps = rep.int(1e-3, length(coefficients)), # A vector of step sizes for the finite-difference approximation to the gradient, on par/parscale scale. Defaults to 1e-3.
                   maxit = 100L, abstol = -Inf,
                   reltol = sqrt(.Machine$double.eps), # Relative convergence tollerance
                   alpha = 1.0, beta = 0.5,
                   gamma = 2.0, REPORT = 10, type = 1, lmm = 5, factr = 1e7,
                   pgtol = 0, tmax = 10, temp = 10.0)
  
  opt_ctrl[names(optim.control)] <- optim.control

  C_coef <- constant_gmm$par

  constants <- lapply(1:(nrow(input_degree) - 1), FUN = function(i) { # input_degree has as many rows as inputs
    new_in_deg <- input_degree # Input degree has first rows for fixed than flex inputs. Here: first row fixed and second flex

    # Take one off from every degree if > 0 and leave at 0 otherwise
    new_in_deg[i, ] <- ifelse(new_in_deg[i, ] > 0,
                              new_in_deg[i, ] - 1,
                              new_in_deg[i, ])
    
    # Only take these where flex input == 0 (last row)
    new_C_deg <- new_in_deg[, new_in_deg[nrow(input_degree), ] == 0]

    C_match <- apply(new_C_deg, MARGIN = 2, FUN = match_gnr, degree_vec =
                       input_degree) # This thing goes through all columns of the derivative of input_degrees and compares it to input_degrees columns.
    
    deriv_C <- all_input[, C_match]
    deriv_C[is.na(deriv_C)] <- 1

    C <- deriv_C %*%
      t(t(input_degree[i, input_degree[nrow(input_degree), ] == 0]) * C_coef)
  })
  
  elas_noC <- lapply(1:(nrow(orig_input_degree) - 1), FUN = function(i) {
    new_in_deg <- orig_input_degree
    new_in_deg[i, ] <- ifelse(new_in_deg[i, ] > 0,
                              new_in_deg[i, ] - 1,
                              new_in_deg[i, ])

    new_in_deg[nrow(new_in_deg), ] <- new_in_deg[nrow(new_in_deg), ] + 1

    in_match <- apply(new_in_deg, MARGIN = 2, FUN = match_gnr, degree_vec =
                        orig_input_degree)

    deriv_input <- orig_input[, in_match]
    deriv_input[is.na(deriv_input)] <- 0
    elas <- deriv_input %*% t(t(orig_input_degree[i, ]) * (object$arg$D_coef))
  })
  print("elas_noC")
  print(elas_noC)
  elas = lapply(1:length(elas_noC), FUN = function(x) {
    elas_noC[[x]] + constants[[x]]
  })

  logomega <- big_Y - (as.matrix(pred) %*% C_coef)
  omega <- exp(logomega)
  productivity <- as.matrix(exp(logomega + errors))

  elasticities <- do.call(cbind, elas)
  colnames(elasticities) <- object$arg$fixed_names
  colnames(productivity) <- "productivity"
  ss_return <- list("fixed_elas" = elasticities,
                    "productivity" = productivity,
                    "optim_info" = constant_gmm,
                    "starting_values" = coefficients,
                    "control" = list("degree_w" = degree_w,
                                     "degree_tau" = degree_tau,
                                     "method" = method,
                                     "optim_control" = opt_ctrl))
  class(ss_return) <- "gnriv"
  return(constants)
}

constant_moments_test <- function(C_kl, data, big_Y_base, big_Y_lag, lag_data,
                             degree) {
  w <- big_Y_base - (data %*% C_kl)
  w_1 <- big_Y_lag - (lag_data %*% C_kl)
  
  if (degree < 2) {
    markov <- w_1
  } else {
    poly <- sapply(2:degree, FUN = function(i) {
      `^`(w_1, i)
    })
    
    markov <- cbind(w_1, poly)
  }

  reg <- stats::lm(w ~ markov)
  csi <- w - reg$fitted.values

  moments <- apply(data, MARGIN = 2, FUN = function(i) {
    sum(i * csi) / length(i)
  })
  
  obj <- t(moments) %*% moments

  # print(markov[1:5])
  # print(length(big_Y_lag))
  # print("taylor_fixed[1]")
  # print(data[1,])  
  # print("lag_data[1]")
  # print(lag_data[1,])
  # print("C_kl")
  # print(C_kl)
  # print("w[1]")
  # print(w[1])
  # print("w_1[1]")
  # print(w_1[1])
  # print(moments)
  # print(obj)
  # print("END")
  # return(obj)
}

constant_moments <- function(C_kl, data, big_Y_base, big_Y_lag, lag_data,
                             degree) {
  w <- big_Y_base - (data %*% C_kl)
  w_1 <- big_Y_lag - (lag_data %*% C_kl)
  
  if (degree < 2) {
    markov <- w_1
  } else {
    poly <- sapply(2:degree, FUN = function(i) {
      `^`(w_1, i)
    })
    
    markov <- cbind(w_1, poly)
  }

  reg <- stats::lm(w ~ markov)
  csi <- w - reg$fitted.values

  moments <- apply(data, MARGIN = 2, FUN = function(i) {
    sum(i * csi) / length(i)
  })

  
  obj <- t(moments) %*% moments

  return(obj)
}

match_gnr <- function(i, degree_vec) {
  match_test <- apply(i == degree_vec, MARGIN = 2, FUN = prod) # For whatever fucking reason, you can use the prod fnc in R to compare vectors.
  col_index <- ifelse(length(which(match_test == 1)) != 0,
                      which(match_test == 1), NA)
}


