model_select <-
function (div, ages, me1 = NULL, me2 = NULL, GRAD =NULL, cats=NULL, breakpoint=NULL, domain=NULL,   
    models = c("BM_null", "OU_null", "BM_linear", "OU_linear", "OU_linear_sig", 
    "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_cat"), 
    starting = NULL, absolute=TRUE, parallel=FALSE, cores=NULL) {
  
  if(length(div)!=length(ages)) {
    stop("Hold up: The vectors for species pair trait differences and species pair age are of unequal length")
  }
  if(is.null(starting)==FALSE) {
      if(length(models)>1 & is.list(starting)==FALSE) {
        stop("Hold up: When assessing more than one model and providing your own starting parameters, 'starting' must \n   be a list in which each element contains starting parameters for one of the models you are testing.")
      }
  }
  if(all(models %in% c("BM_null", "OU_null", "BM_linear", "OU_linear", "OU_linear_sig", 
    "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_cat")) == FALSE) {
    stop("Spell check: you've entered at least one model that doesn't match the models accepted by this function")
  }
  if(sum(is.null(me1), is.null(me2)) == 1) stop("You've supplied measurement error 
    for just one species in the pairs; please provide 2 or 0 vectors for measurement error")

  # create shared output matrix to summarize results from likelihood searches
  RESULTS_SUMMARY <- matrix(NA, length(models), 19)
  colnames(RESULTS_SUMMARY) <- c("model", "logLik", "k (no. params)", "AIC", 
        "AICc", "delta_AICc", "Akaike Weight", "alpha", "alpha_slope", "sig2", "sig2_slope", "psi_int", "psi_slope", 
        "psi2_int", "psi2_slope", "wait_time", "convergence", "psi2", "psi3") 
  rownames(RESULTS_SUMMARY) <- models
  N <- length(models)
  if (length(starting) == 1) {
    if (is.null(starting[1])) starting[1:N] = NULL
  }

  # select starting parameters and run the estimation for models chosen by user
  if ("BM_null" %in% models) { 
    j <- which(models == "BM_null")  
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "BM_null", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "BM_null", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1) 
    RESULTS_SUMMARY[j, 10] <- res$par[1] # beta
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }
  if ("BM_linear" %in% models) {
    j <- which(models == "BM_linear")  
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "BM_linear", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "BM_linear", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 10] <- res$par[2] # beta int
    RESULTS_SUMMARY[j, 11] <- res$par[1] # beta slope
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }
  if ("OU_null" %in% models) {
    j <- which(models == "OU_null")  
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "OU_null", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "OU_null", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 8] <- res$par[1] # alpha
    RESULTS_SUMMARY[j, 10] <- res$par[2] # beta 
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }
  if ("OU_linear" %in% models) {
    j <- which(models == "OU_linear")  
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "OU_linear", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "OU_linear", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 8] <- res$par[1] # alpha int
    RESULTS_SUMMARY[j, 9] <- res$par[2] # alpha slope
    RESULTS_SUMMARY[j, 10] <- res$par[3] # beta 
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }
  if ("OU_linear_sig" %in% models) {
    j <- which(models == "OU_linear_sig")  
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "OU_linear_sig", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "OU_linear_sig", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 8] <- res$par[1] # alpha
    RESULTS_SUMMARY[j, 11] <- res$par[2] # beta slope
    RESULTS_SUMMARY[j, 10] <- res$par[3] # beta int
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }
  if ("DA_null" %in% models) {
    j <- which(models == "DA_null")  
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_null", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_null", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 8] <- res$par[1] # alpha constant
    RESULTS_SUMMARY[j, 10] <- res$par[2] # beta constant
    RESULTS_SUMMARY[j, 12] <- res$par[3] # psi constant
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }
  if ("DA_linear" %in% models) {
    if(is.null(GRAD)) {
      stop("Hold up! You've chosen to assess a model that requires the gradient position of each pair as input")
    }
    j <- which(models == "DA_linear")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_linear", p_starting = starting[[j]], 
        div = div, ages= ages, me1 = me1, me2= me2, GRAD = GRAD, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_linear", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 8] <- res$par[1] # alpha constant
    RESULTS_SUMMARY[j, 10] <- res$par[2] # beta constant
    RESULTS_SUMMARY[j, 12] <- res$par[4] # psi_int
    RESULTS_SUMMARY[j, 13] <- res$par[3] # psi_slope
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }
  if ("DA_cat" %in% models) {
    if(all(cats <= 2)==FALSE) {
      stop("Hold up! This function can only be used for up to three categorical variables, which must 
        be coded \n with cat values of 0, 1, and 2, respectively")
    }
    if(is.null(cats)) {
      stop("Hang on, you have to provide the cats vector")
    }
    j <- which(models == "DA_cat")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_cat", p_starting = starting[[j]], 
        div = div, ages= ages, me1 = me1, me2= me2, GRAD = GRAD, cats=cats, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_cat", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, cats=cats, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 8] <- res$par[1] # alpha constant
    RESULTS_SUMMARY[j, 10] <- res$par[2] # beta constant
    RESULTS_SUMMARY[j, 12] <- res$par[3] # psi1
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
    RESULTS_SUMMARY[j, 18] <- res$par[4] # psi2
    if(length(res$par)==5) RESULTS_SUMMARY[j, 19] <- res$par[5] # psi3
  }
  if ("DA_wt" %in% models) {
    j <- which(models == "DA_wt")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_wt", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_wt", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par)# no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 8] <- res$par[1] # alpha constant
    RESULTS_SUMMARY[j, 10] <- res$par[2] # beta constant
    RESULTS_SUMMARY[j, 12] <- res$par[3] # psi1 constant
    RESULTS_SUMMARY[j, 14] <- res$par[4] # psi2 constant
    RESULTS_SUMMARY[j, 16] <- res$par[5] # wait time 
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }
  if ("DA_bp" %in% models) {
    j <- which(models == "DA_bp")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_bp", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, breakpoint=breakpoint, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_bp", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, breakpoint=breakpoint, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 8] <- res$par[1] # alpha constant
    RESULTS_SUMMARY[j, 10] <- res$par[2] # beta constant
    RESULTS_SUMMARY[j, 12] <- res$par[3] # psi1_int 
    RESULTS_SUMMARY[j, 14] <- res$par[4] # psi2_int
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }
  if ("DA_wt_linear" %in% models) {
    j <- which(models == "DA_wt_linear")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_wt_linear", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_wt_linear", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2]) 
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 8] <- res$par[1] # alpha constant
    RESULTS_SUMMARY[j, 10] <- res$par[2] # beta constant
    RESULTS_SUMMARY[j, 12] <- res$par[4] # psi1_int 
    RESULTS_SUMMARY[j, 13] <- res$par[5] # psi1_slope
    RESULTS_SUMMARY[j, 14] <- res$par[6] # psi2_int
    RESULTS_SUMMARY[j, 15] <- res$par[7] # psi2_slope
    RESULTS_SUMMARY[j, 16] <- res$par[3] # wait time
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }
  if ("DA_bp_linear" %in% models) {
    j <- which(models == "DA_bp_linear")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_bp_linear", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, breakpoint=breakpoint, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_bp_linear", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, breakpoint=breakpoint, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    RESULTS_SUMMARY[j, 2] <- -res$objective # log likelihood
    RESULTS_SUMMARY[j, 3] <- length(res$par) # no. params
    RESULTS_SUMMARY[j, 4] <- 2 * as.numeric(RESULTS_SUMMARY[j,3]) - 2 * as.numeric(RESULTS_SUMMARY[j, 2])
    RESULTS_SUMMARY[j, 5] <- as.numeric(RESULTS_SUMMARY[j,4]) + (2 * as.numeric(RESULTS_SUMMARY[j, 3]) * (as.numeric(RESULTS_SUMMARY[j,3]) + 1))/(length(div) - as.numeric(RESULTS_SUMMARY[j, 3]) - 1)
    RESULTS_SUMMARY[j, 8] <- res$par[1] # alpha constant
    RESULTS_SUMMARY[j, 10] <- res$par[2] # beta constant
    RESULTS_SUMMARY[j, 12] <- res$par[4] # psi1_int 
    RESULTS_SUMMARY[j, 13] <- res$par[5] # psi1_slope
    RESULTS_SUMMARY[j, 14] <- res$par[6] # psi2_int
    RESULTS_SUMMARY[j, 15] <- res$par[7] # psi2_slope
    RESULTS_SUMMARY[j, 17] <- res$convergence # termination code from nlminb
  }

  # summarize results in output matrix
  RESULTS_SUMMARY[which(RESULTS_SUMMARY[,17]>0),5]=NA # converts AICc to NA for models in which convergence wasn't succesful
  RESULTS_SUMMARY[,6] = RESULTS_SUMMARY[,5] - min(RESULTS_SUMMARY[which(is.na(RESULTS_SUMMARY[,5])==FALSE),5]) # delta AICc calculated from only those models with succesful convergence
  RESULTS_SUMMARY[,7] = exp(-0.5*RESULTS_SUMMARY[,6])/sum(exp(-0.5*RESULTS_SUMMARY[,6]),na.rm=TRUE) # Akaike weights calc'd from deltAICc
  RESULTS_SUMMARY <- RESULTS_SUMMARY[, -1, drop = FALSE]
  RESULTS_SUMMARY <- t(RESULTS_SUMMARY)
  RESULTS_SUMMARY=apply(X=RESULTS_SUMMARY, MARGIN=2, FUN=round, digits=5)
  return(RESULTS_SUMMARY)
}
