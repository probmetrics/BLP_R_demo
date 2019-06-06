##
## Functions for estimating market share in BLP model
##
## Shu Xu, <May, 2019>
##

share_per_mkt <- function(delta, alpha, sigma, p, X, y, vshk) {
    
    # 
    # ** Estimating market share in a single market **
    # 
    # Inputs: 
    #   delta,  -- J vector
    #   X,      -- J x K matrix
    #   y,      -- S vector of income y_i
    #   p,      -- J vector of price p_j 
    #   vshk,   -- K x S matrix
    #   alpha,  -- scalar, coef of price/income
    #   sigma,  -- K vector
    #
    # Output:
    #   J x 1 vector of estimated market share
    #
    
    require(matrixStats)
    
    # ubar <- delta - alpha * outer(p, y, "/") + X %*% sigma %*% vshk     
    ubar <- delta - alpha * outer(p, y, "/") + X %*% (sigma * vshk) #<-- J by S matrix
    ubars <- cbind(0, t(ubar)) #<-- S x (J + 1) matrix
    
    prob <- exp(ubars - rowLogSumExps(ubars))[, -1]
    share <- colMeans(prob)
    return(share)
}


mkt_share_est <- function(mkt_idx, delta, alpha, sigma, price, X, yshk, vshk){
    
    # 
    # ** Estimating market share for all markets **
    # 
    # Inputs: 
    #   mkt_idx  -- index indicating diff. markets
    #   delta,   -- J x M vector
    #   X,       -- JM x K matrix
    #   yshk     -- M x S matrix of income y_i
    #   p,       -- JM vector of price p_j 
    #   vshk,    -- (K x S x M) 3D array
    #   alpha,   -- scalar, coef of price/income
    #   sigma,   -- K vector
    #
    # Output:
    #   JM x 1 vector of estimated market shares
    #
    
    mid <- unique(mkt_idx)
    
    mshr_fun <- function(i) {
        midx <- mkt_idx == mid[i]
        mshr <- share_per_mkt(delta[midx], alpha, sigma, price[midx], X[midx, ], yshk[i, ], vshk[, , i])
        return(mshr)
    }

    share_list <- sapply(1:length(mid), mshr_fun)       
    mkt_share <- unlist(share_list)
    
    return(mkt_share)
}
