
##
## Prepare BLP Data
##

setwd("E:/NutStore/My Lecture Scripts/PhD_Frontiers_in_Applied_Micro/Code_Examples/BLP_funs")

# library(hdm)
# data("BLP")
# BLP_data <- BLP$BLP
# BLP_data <- transform(BLP_data, price = price + 11.761)

BLP_data <- read.csv("BLP_data.csv", header = TRUE) 

# --- identical firm index maxtrix within a market ---
mkt_idx <- BLP_data$cdid
fid_list <- sapply(unique(mkt_idx), function(m){
    fid <- BLP_data$firm.id[mkt_idx == m]
    outer(fid, fid, function(x, y) as.numeric(x == y))
})


# --- product characteristic matrix ---
X <- as.matrix(cbind(cons = 1, BLP_data[, c("hpwt", "air", "mpd", "space")]))
W <- as.matrix(cbind(cons = 1, ln.hpwt = log(BLP_data$hpwt), air = BLP_data$air, 
                     ln.mpg = log(BLP_data$mpg),ln.space = log(BLP_data$space), 
                     trend = BLP_data$trend))
price <- BLP_data$price
data_share <- BLP_data$share

##
## Draw simulated v, y
## From IPUMS-CPS average log (income in thousands of 1983 $) = 3.082,
## with little variation across years, sd = 0.840
##

# --- random draws ---
library(randtoolbox)
nsim <- 200
nmkt <- length(unique(mkt_idx))
npar <- ncol(X)

eshk <- halton(nsim * nmkt, dim = npar + 1, normal = TRUE)
yshk <- matrix(exp(eshk[, npar + 1]*0.840 + 3.082), nrow = nmkt)
vshk <- array(eshk[, 1:npar], dim = c(npar, nsim, nmkt))

# --- setup parameters ---
delta_tmp <- rnorm(nrow(BLP_data), mean = 1, sd = 0.1)
delta_init <- with(BLP_data, log(share) - log(outshr))

alpha <- 43.501
sigma <- c(3.612, 4.628, 1.818, 1.050, 2.056)
beta <- c(-7.061, 2.883, 1.521, -0.122, 3.460)
gamma <- c(0.952, 0.477, 0.619, -.415, -.049, .019)

##
## Testing functions
##

source("mkt_share_est.R")

# --- estimating market share ---
pred_share <- mkt_share_est(mkt_idx, delta_init, alpha, sigma, price, X, yshk, vshk)

# --- contraction mapping ---
delta_sem <- delta_fpt_sem(mkt_idx, delta_init, alpha, sigma, data_share, price, X, yshk, vshk)
delta_cmp <- delta_fpt_cmp(mkt_idx, delta_init, alpha, sigma, data_share, price, X, yshk, vshk)

# --- log marginal cost ---
lnmc <- marg_cost_est(mkt_idx, delta_cmp, alpha, sigma, data_share, price, X, yshk, vshk, fid_list)      

# --- linear parameters ---
# demand side IV 
Zd <- hdm:::constructIV(BLP_data$firm.id, BLP_data$cdid, BLP_data$id, X)
Zd <- as.matrix(cbind(X, Zd))

# supply side IV
Zs <- hdm:::constructIV(BLP_data$firm.id, BLP_data$cdid, BLP_data$id, W)
Zs <- as.matrix(cbind(W, Zs))

# --- first step GMM ---
lgmm_fit(delta_sem, X, Zd)$coefs
lgmm_fit(lnmc, W, Zs)$coefs

yds <- c(delta_sem, lnmc)
Xds <- Matrix::bdiag(X, W)
colnames(Xds) <- c(colnames(X), colnames(W))
Zds <- Matrix::bdiag(Zd, Zs)
colnames(Zds) <- c(colnames(Zd), colnames(Zs))
lgmm_fit(yds, Xds, Zds)$coefs


                      