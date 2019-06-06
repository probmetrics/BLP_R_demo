##
## Functions for calculating price elasticity and (log) marginal cost in BLP model
##
## Shu Xu, <May, 2019>
##


dshare_dp <- function(delta, alpha, sigma, p, X, y, vshk, fid_mat) {
    
    # 
    # ** Matrix of market share derivative w.r.t price (dmkt_shr/dp) **
    # 
    # Inputs: 
    #   delta,      -- J vector of fixed point delta!
    #   X,          -- M x K matrix
    #   y           -- S vector of income y_i
    #   p,          -- J vector of price p_j 
    #   vshk,       -- K x S matrix
    #   alpha,      -- scalar, coef of price/income
    #   sigma,      -- K vector
    #   fid_mat     -- J x J matrix, indicating whether two products come from same firm
    #
    # Output:
    #   J x J matrix of price elasticity
    #
    
    require(matrixStats)
    
    ubar <- delta - alpha * outer(p, y, "/") + X %*% (sigma * vshk)
    ubars <- cbind(0, t(ubar)) #<-- S x (J + 1) matrix
    prob <- exp(ubars - rowLogSumExps(ubars))[, -1] #<-- S x J matrix
    
    dmu_dp <- - alpha / y      #<-- S x 1
    share.dmu_dp <- colMeans(prob * dmu_dp)
    
    dshare_full <- diag(share.dmu_dp) - crossprod(prob, prob * dmu_dp) / nrow(prob)
    dshare <- dshare_full * fid_mat

    return(dshare)
}


marg_cost_est <- function(mkt_idx, delta, alpha, sigma, data_share, 
                          price, X, yshk, vshk, fid_list) {
    
    # 
    # ** Vector of (log) marginal cost **
    #
    # Inputs: 
    #   mkt_idx     -- index indicating diff. markets
    #   delta,      -- JM vector of fixed point delta!
    #   X,          -- JM x K matrix
    #   yshk        -- M x S matrix of income y_i
    #   price,      -- JM vector of price p_j 
    #   vshk,       -- K x S x M Array
    #   alpha,      -- scalar, coef of price/income
    #   sigma,      -- K vector of random shock's std
    #   fid_list    -- M-length list, each element is a  J x J fid_mat matrix
    #
    # Output:
    #   JM x 1 vector of (log) marginal cost
    #
    
    mid <- unique(mkt_idx)
    
    mcfun <- function(i) {
        midx <- mkt_idx == mid[i]
        dsdp <- dshare_dp(delta[midx], alpha, sigma, price[midx], X[midx, ], 
                          yshk[i, ], vshk[, , i], fid_list[[i]])
        
        lnmc <- log(price[midx] - solve(dsdp) %*% data_share[midx])
        return(lnmc)
    }
    
    lnmc_list <- sapply(1:length(mid), mcfun)
    return(unlist(lnmc_list))
}
