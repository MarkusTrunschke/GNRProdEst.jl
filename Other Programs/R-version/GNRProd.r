# Get newest version of GNRProd (R pakcage from DavidJ.Jin<jindavid@sas.upenn.edu>)
devtools::install_github("davidjin0/gnrprod")
library(gnrprod)

######################### Test GNRProd with Columbian data (as David shows on Github) #########################
# load Colombian plant-level data
# Columbia_311 <- colombian
# write.csv(Columbia_311, "C:/Users/marku/Documents/GNRProdEst/Other Programs/R-version/Columbia_311.csv")
Columbia_311 <- read.csv("C:/Users/marku/Documents/GNRProdEst/Other Programs/R-version/Columbia_311.csv")

# estimate production function parameters and productivity
gnr_fit <- gnrprod(output = "RGO",
                   fixed = c("L", "K"),
                   flex = "RI",
                   share = "share",
                   id = "id",
                   time = "year",
                   data = Columbia_311)
summary(gnr_fit) # It works

gnr_FS <- gnrflex(output = "RGO",
                   fixed = c("L", "K"),
                   flex = "RI",
                   share = "share",
                   id = "id",
                   time = "year",
                   data = Columbia_311,
                   control = list(degree = 2, maxit = 2000))


int_G_I2 = gnr_FS$integ_G_I
flex_elas_list = gnr_FS$elas
flex_elas = flex_elas_list$flex_elas
fs_coefs = flex_elas_list$coefficients

######################### Now use GNR Monte Carlo data (Cobb-Douglas version to stay consistent over examples) #########################
GNR_MC_data = read.fwf(file = "C:/Users/marku/Documents/GNRProdEst/Other Programs/GNR/Table_1/Monte_Carlo_Data/cd_data.out", widths = c(11,11,11,11,11,11))
# Rename variables (deducted from "our_estimator.do" file)
names(GNR_MC_data)[1] <- "ID"
names(GNR_MC_data)[2] <- "Year"
names(GNR_MC_data)[3] <- "h"
names(GNR_MC_data)[4] <- "K"
names(GNR_MC_data)[5] <- "M"
names(GNR_MC_data)[6] <- "Y"

write.csv(GNR_MC_data, "C:/Users/marku/Documents/GNRProdEst/Other Programs/R-version/cd_data.csv")
GNR_MC_data <- read.csv("C:/Users/marku/Documents/GNRProdEst/Other Programs/GNR/Cleaned version of Table_1/cd_data_500.csv")

# Calculate share variable
GNR_MC_data$m_y_share <- GNR_MC_data$i_level / GNR_MC_data$yg_level # nolint log(share of intermediate input to revenue variable)
GNR_MC_data$ln_m_y_share <- log(GNR_MC_data$i_level / GNR_MC_data$yg_level) # nolint log(share of intermediate input to revenue variable)

test = m_y_share <- GNR_MC_data$M / GNR_MC_data$Y

# GNR use batches of 500 firms in their MC. Do the same
GNR_MC_sample <- GNR_MC_data[GNR_MC_data$id < 500, ]

# Run GNR estimator
gnr_est <- gnrprod(data = GNR_MC_sample,
                   output = "y",
                   fixed = "k",
                   flex = "m",
                   id = "ID",
                   time = "Year",
                   share = "ln_m_y_share",
                   ss_control = list(method = "BFGS", trace = 1, maxit = 2000)) # nolint
summary(gnr_est) # Converges and gives decent results

GNR_MC_sample$big_Y = gnr_FS$arg$big_Y
GNR_MC_sample$errors <- gnr_FS$elas$residuals


write.csv(GNR_MC_sample, "C:/Users/marku/Documents/GNRProdEst/Other Programs/GNR/Cleaned version of Table_1/cd_data_500_w_fs.csv")

gnr_FS <- gnrflex(output = "yg",
                   fixed = "k",
                   flex = "i",
                   id = "id",
                   time = "time",
                   share = "si",
                   data = GNR_MC_sample,
                   control = list(degree = 2, maxit = 2000))

gnr_SES <- gnriv(object = gnr_FS)

fs_elas_list = gnr_FS$elas
fs_arg_list = gnr_FS$arg

big_Y = fs_arg_list$big_Y
epsilon = fs_elas_list$residuals

flex_elas = fs_elas_list$flex_elas
fs_coefs = fs_elas_list$coefficients
######################### Now try my generated data #########################
test_data = read.csv(file = "C:/Users/marku/Documents/GNRProdEst/Data/test_data_CD_K_L_M.csv") # nolint

test_data$m_y_share <- log(test_data$M / test_data$Y) # nolint log(share of intermediate input to revenue variable)

test_sample <- test_data[test_data$ID < 2000, ]

# Run GNR estimator
gnr_est <- gnrprod(data = test_sample,
                   output = "y",
                   fixed = c("k", "l"),
                   flex = "m",
                   id = "ID",
                   time = "year",
                   share = "m_y_share",
                   ss_control = list(method = "BFGS", trace = 1, maxit = 20000)) # nolint
summary(gnr_est) # Converges but gives results that are far from the true parameters