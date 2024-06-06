# Get newest version of GNRProd (R pakcage from DavidJ.Jin<jindavid@sas.upenn.edu>)
devtools::install_github("davidjin0/gnrprod")
library(gnrprod)

######################### Test GNRProd with Columbian data (as David shows on Github) #########################
# load Colombian plant-level data
data <- colombian

# write.csv(data, "/Users/markus_trunschke/Documents/GNRProdEst/Other Programs/R-version/colombian.csv")
data <- read.csv(file = "/Users/markus_trunschke/Documents/GNRProdEst/Other Programs/R-version/colombian.csv")
# estimate production function parameters and productivity
gnr_fit <- gnrprod(output = "RGO",
                   fixed = c("L", "K"),
                   flex = "RI",
                   share = "share",
                   id = "id",
                   time = "year",
                   data = data)
summary(gnr_fit) # It works

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

# Calculate logs and share variable
GNR_MC_data$y <- log(GNR_MC_data$Y)
GNR_MC_data$k <- log(GNR_MC_data$K)
GNR_MC_data$m <- log(GNR_MC_data$M)
GNR_MC_data$m_y_share <- log(GNR_MC_data$M / GNR_MC_data$Y) # nolint log(share of intermediate input to revenue variable)

# GNR use batches of 500 firms in their MC. Do the same
GNR_MC_sample <- GNR_MC_data[GNR_MC_data$ID < 500, ]

# Run GNR estimator
gnr_est <- gnrprod(data = GNR_MC_sample,
                   output = "y",
                   fixed = "k",
                   flex = "m",
                   id = "ID",
                   time = "Year",
                   share = "m_y_share",
                   ss_control = list(method = "BFGS", trace = 1, maxit = 2000)) # nolint
summary(gnr_est) # Converges and gives decent results


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