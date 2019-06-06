
##
## Functions for doing contraction mapping in BLP model
##
##  Shu Xu, <May, 2019>
##

delta_fpt_sem <- function(mkt_idx, delta_init, alpha, sigma, data_share, price,  
                          X, yshk, vshk, tol = 1E-12, trace = TRUE, ...) {

    # 
    # ** Contraction mapping to find delta using SQUAREM **
    # 
    # Note: delta(s) found by contraction mapping is done for all markets at once.
    #
    # Inputs: 
    #   mkt_idx     -- index indicating diff. markets
    #   delta_init, -- JM vector
    #   X,          -- JM x K matrix
    #   yshk        -- M x S matrix of income y_i
    #   p,          -- JM vector of price p_j 
    #   vshk,       -- K x S x M array
    #   alpha,      -- scalar, coef of price/income
    #   sigma,      -- K vector
    #
    # Output:
    #   JM x 1 vector of fixed point delts(s)
    #
    
    require(SQUAREM)
    
    # -- define the fixed point function --
    fixptfn <- function(delta_old) {
        pred_share <- mkt_share_est(mkt_idx, delta_old, alpha, sigma, price, X, yshk, vshk)
        pred_share[pred_share <= 0] <- min(pred_share[pred_share > 0])
        
        delta_new <- delta_old + log(data_share) - log(pred_share)
        delta_old <- delta_new
        return(delta_new)
    }
    
    # -- do the contraction mapping --
    out <- squarem(delta_init, fixptfn, control = list(tol = tol, trace = trace, ...))
    return(out$par)
}

delta_fpt_cmp <- function(mkt_idx, delta_init, alpha, sigma, data_share, price,  
                          X, yshk, vshk, print.skip = 5, niters = 1000, tol = 1E-12){
    
    # 
    # ** Contraction mapping to find delta in traditional way **
    # 
    # Note: delta(s) found by contraction mapping is done for all markets at once.
    #
    # Inputs: 
    #   mkt_idx     -- index indicating diff. markets
    #   delta_init, -- JM vector
    #   X,          -- JM x K matrix
    #   yshk        -- M x S matrix of income y_i
    #   p,          -- JM vector of price p_j 
    #   vshk,       -- K x S x M array
    #   alpha,      -- scalar, coef of price/income
    #   sigma,      -- K vector
    #
    # Output:
    #   JM x 1 vector of fixed point delts(s)
    #
    
    err <- 10
    iter <- 0
    delta_old <- delta_init
    while ((err > tol) & (iter < niters)) {
        
        pred_share <- mkt_share_est(mkt_idx, delta_old, alpha, sigma, price, X, yshk, vshk)
        pred_share[pred_share <= 0] <- min(pred_share[pred_share > 0])
        delta_new <- delta_old + log(data_share) - log(pred_share)
        
        err <- max(abs(delta_old - delta_new))
        delta_old <- delta_new
        iter <- iter + 1
        
        if (iter %% print.skip == 1 | err <= tol) {
            cat("the", iter, "th iterations,", "rel.tol =", err, "\n")
        } 
        
    }
    
    if(iter >= niters) { 
        stop("Maximum iters reached. The number of Contraction Mapping iteration needs to be increased.\n") 
    }
    
    return(delta_new)
}
