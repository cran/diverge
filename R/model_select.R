model_select <-
function (div, ages, me1 = NULL, me2 = NULL, GRAD =NULL, cats=NULL, breakpoint=NULL, domain=NULL,   
    models, starting = NULL, absolute=TRUE, parallel=FALSE, cores=NULL) {
  
# warnings
  if(length(div)!=length(ages)) {
    stop("Hold up: The vectors for species pair trait differences and species pair age are of unequal length")
  }
  if(any(is.na(div)) | any(is.na(ages))) {
    stop("It looks like you've got missing values for pair age or trait differences - please remove the missing indices from all input vectors")
  }
  if(any(is.na(GRAD))) {
    stop("It looks like you've got missing values for the gradient variable - please remove the missing value(s) and remove the corresponding values for pair age and trait divergence (and any other input vector)")
  }
  if(is.null(starting)==FALSE) {
      if(length(models)>1 & is.list(starting)==FALSE) {
        stop("Hold up: When assessing more than one model and providing your own starting parameters, 'starting' must \n   be a list in which each element contains starting parameters for one of the models you are testing.")
      }
  }
  if(all(models %in% c("BM_null", "BM_linear", "BM_cat", "OU_null", "OU_linear", "OU_linear_sig", 
    "OU_cat", "DA_null", "DA_linear", "DA_cat", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_OU", "DA_OU_linear", 
    "DA_OU_cat", "DA_BM", "DA_BM_linear", "DA_BM_cat", "OU_BM", "OU_BM_linear", "OU_BM_cat")) == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }
  if(sum(is.null(me1), is.null(me2)) == 1) stop("You've supplied measurement error 
    for just one species in the pairs; please provide either 2 or 0 vectors for measurement error")

# create output matrix
  res.summ = matrix(NA, nrow=10, ncol=length(models))
  rownames(res.summ) <- c("logLik", "k (no. params)", "AIC", 
        "AICc", "delta_AICc", "Akaike Weight", "convergence", "alpha", "sig2", "psi") 
  colnames(res.summ) <- models
  if(is.null(GRAD)==FALSE) {
    res.summ = rbind(res.summ, "alpha_slope"=rep(NA,length(models)), 
      "sig2_slope"=rep(NA,length(models)), "psi_slope"=rep(NA,length(models)))
  }
  if("BM_cat" %in% models | "OU_cat" %in% models | "DA_cat" %in% models | "DA_wt" %in% models 
    | "DA_wt_linear" %in% models | "DA_bp" %in% models | "DA_bp_linear" %in% models) {
    parmat = matrix(NA, ncol=length(models), nrow=(length(unique(cats))-1)*3)
    rownames(parmat) = c(paste("alpha", 2:length(unique(cats)), sep="_"), 
      paste("sig2", 2:length(unique(cats)), sep="_"), paste("psi", 2:length(unique(cats)), sep="_"))
    res.summ = rbind(res.summ, parmat)
  }
  if("DA_wt" %in% models | "DA_wt_linear" %in% models) {
    res.summ = rbind(res.summ, "wait_time"=rep(NA,length(models)))
  }
  if("DA_wt_linear" %in% models | "DA_bp_linear" %in% models) {
    res.summ = rbind(res.summ, "psi2_slope"=rep(NA, length(models)))
  }
  if("DA_OU" %in% models | "DA_BM" %in% models | "OU_BM" %in% models) {
    res.summ=rbind(res.summ,"prop"=rep(NA, length(models)))
  }
  if("DA_OU_linear" %in% models | "DA_BM_linear" %in% models | "OU_BM_linear" %in% models) {
    res.summ=rbind(res.summ,"prop"=rep(NA, length(models)),"prop_slope"=rep(NA, length(models)))
  }
  if("DA_OU_cat" %in% models | "DA_BM_cat" %in% models | "OU_BM_cat" %in% models) {
    if("DA_OU" %in% models == FALSE & "DA_BM" %in% models == FALSE & "OU_BM" %in% models == FALSE) {
      stop("You have to include DA_OU, DA_BM, or OU_BM model in addition to the categorical model(s) you've chosen")
    }
    parmat2 = matrix(NA, ncol=length(models), nrow=(length(unique(cats))-1))
    rownames(parmat2) = c(paste("prop", 2:length(unique(cats)), sep="_"))
    res.summ=rbind(res.summ, parmat2)
  }

# select starting parameters and run the ML estimation for models chosen by user
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1) 
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[9, j] <- res$par[1] # sig2
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[9, j] <- res$par[2] # sig2 int
    res.summ[12, j] <- res$par[1] # sig2 slope
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3,j]) + (2 * as.numeric(res.summ[2,j]) * (as.numeric(res.summ[2,j]) + 1))/(length(div) - as.numeric(res.summ[2,j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha
    res.summ[9, j] <- res$par[2] # sig2
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
    res.summ[1,j] <- -res$objective # log likelihood
    res.summ[2,j] <- length(res$par) # no. params
    res.summ[3,j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4,j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7,j] <- res$convergence # termination code from nlminb
    res.summ[8,j] <- res$par[1] # alpha int
    res.summ[9,j] <- res$par[3] # sig2
    res.summ[11, j] <- res$par[2] # alpha slope
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha
    res.summ[9, j] <- res$par[3] # sig2 int
    res.summ[12, j] <- res$par[2] # sig2 slope
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[3] # psi constant
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2 ,j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # beta constant
    res.summ[10, j] <- res$par[4] # psi_int
    res.summ[13, j] <- res$par[3] # psi_slope
  }
  if ("DA_cat" %in% models) {
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[3] # psi1
    for(i in 2:length(unique(cats))) {
      res.summ[paste("psi", i, sep="_"), j] = res$par[3+(i-1)]
    }
  }
  if ("BM_cat" %in% models) {
    if(is.null(cats)) {
      stop("Hang on, you have to provide the cats vector")
    }
    j <- which(models == "BM_cat")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "BM_cat", p_starting = starting[[j]], 
        div = div, ages= ages, me1 = me1, me2= me2, GRAD = GRAD, cats=cats, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "BM_cat", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, cats=cats, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[9, j] <- res$par[1] # sig2_1 constant
    for(i in 2:length(unique(cats))) {
      res.summ[paste("sig2", i, sep="_"), j] = res$par[1+(i-1)]
    }
  }
  if ("OU_cat" %in% models) {
    if(is.null(cats)) {
      stop("Hang on, you have to provide the cats vector")
    }
    j <- which(models == "OU_cat")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "OU_cat", p_starting = starting[[j]], 
        div = div, ages= ages, me1 = me1, me2= me2, GRAD = GRAD, cats=cats, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "OU_cat", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, cats=cats, domain=domain, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha1
    res.summ[9, j] <- res$par[2] # sig2 constant
    for(i in 2:length(unique(cats))) {
      res.summ[paste("alpha", i, sep="_"), j] = res$par[2+(i-1)]
    }
  }
  if ("DA_OU" %in% models) {
    j <- which(models == "DA_OU")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_OU", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_OU", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par)# no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[3] # psi constant
    res.summ["prop", j] <- res$par[4] # proportion under DA
  }
  if ("DA_OU_linear" %in% models) {
    j <- which(models == "DA_OU_linear")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_OU_linear", p_starting = starting[[j]], domain=domain,
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_OU_linear", p_starting = NULL, div = div, domain=domain,
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[3] # psi constant
    res.summ["prop_slope", j] <- res$par[4] # slope of proportion DA
    res.summ["prop", j] <- res$par[5] # intercept of proportion DA
  }
  if ("DA_OU_cat" %in% models) {
    if(is.null(cats)) {
      stop("Hang on, you have to provide the cats vector")
    }
    j <- which(models == "DA_OU_cat")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_OU_cat", p_starting = starting[[j]], 
        div = div, ages= ages, me1 = me1, me2= me2, cats=cats, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_OU_cat", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, cats=cats, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[3] # psi constant
    res.summ["prop", j] <- res$par[4] # proportion under DA in cat 1
    for(i in 2:length(unique(cats))) {
      res.summ[paste("prop", i, sep="_"), j] = res$par[4+(i-1)]
    }
  }
  if ("DA_BM" %in% models) {
    j <- which(models == "DA_BM")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_BM", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_BM", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par)# no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[3] # psi constant
    res.summ["prop", j] <- res$par[4] # proportion under DA
  }
  if ("DA_BM_linear" %in% models) {
    j <- which(models == "DA_BM_linear")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_BM_linear", p_starting = starting[[j]], domain=domain,
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_BM_linear", p_starting = NULL, div = div, domain=domain,
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[3] # psi constant
    res.summ["prop_slope", j] <- res$par[4] # slope of proportion DA
    res.summ["prop", j] <- res$par[5] # intercept of proportion DA
  }
  if ("DA_BM_cat" %in% models) {
    if(is.null(cats)) {
      stop("Hang on, you have to provide the cats vector")
    }
    j <- which(models == "DA_BM_cat")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "DA_BM_cat", p_starting = starting[[j]], 
        div = div, ages= ages, me1 = me1, me2= me2, cats=cats, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "DA_BM_cat", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, cats=cats, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[3] # psi constant
    res.summ["prop", j] <- res$par[4] # proportion under DA in cat 1
    for(i in 2:length(unique(cats))) {
      res.summ[paste("prop", i, sep="_"), j] = res$par[4+(i-1)]
    }
  }
  if ("OU_BM" %in% models) {
    j <- which(models == "OU_BM")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "OU_BM", p_starting = starting[[j]], 
        div = div, ages = ages, me1 = me1, me2= me2, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "OU_BM", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par)# no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ["prop", j] <- res$par[3] # proportion under OU
  }
  if ("OU_BM_linear" %in% models) {
    j <- which(models == "OU_BM_linear")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "OU_BM_linear", p_starting = starting[[j]], domain=domain,
        div = div, ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "OU_BM_linear", p_starting = NULL, div = div, domain=domain,
        ages = ages, me1 = me1, me2= me2, GRAD = GRAD, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ["prop_slope", j] <- res$par[3] # slope of proportion OU
    res.summ["prop", j] <- res$par[4] # intercept of proportion OU
  }
  if ("OU_BM_cat" %in% models) {
    if(is.null(cats)) {
      stop("Hang on, you have to provide the cats vector")
    }
    j <- which(models == "OU_BM_cat")
    if (is.null(starting[[j]][1]) == FALSE) {
      res <- find_mle(model = "OU_BM_cat", p_starting = starting[[j]], 
        div = div, ages= ages, me1 = me1, me2= me2, cats=cats, absolute=absolute, parallel=parallel, cores=cores)
    }
    if (is.null(starting[[j]][1])) {
      res <- find_mle(model = "OU_BM_cat", p_starting = NULL, div = div, 
        ages = ages, me1 = me1, me2= me2, cats=cats, absolute=absolute, parallel=parallel, cores=cores)
    }
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par) # no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ["prop", j] <- res$par[3] # proportion under OU in cat 1
    for(i in 2:length(unique(cats))) {
      res.summ[paste("prop", i, sep="_"), j] = res$par[3+(i-1)]
    }
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par)# no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[3] # psi1 constant
    res.summ["psi_2", j] <- res$par[4] # psi2 constant
    res.summ["wait_time", j] <- res$par[5] # wait time 
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par)# no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[3] # psi1 constant
    res.summ["psi_2", j] <- res$par[4] # psi2 constant
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par)# no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[5] # psi1 constant
    res.summ["psi_slope", j] <- res$par[4] # psi1_slope
    res.summ["psi_2", j] <- res$par[7] # psi2 constant
    res.summ["wait_time", j] <- res$par[3] # wait time 
    res.summ["psi2_slope", j] <- res$par[6] # psi2_slope
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
    res.summ[1, j] <- -res$objective # log likelihood
    res.summ[2, j] <- length(res$par)# no. params
    res.summ[3, j] <- 2 * as.numeric(res.summ[2, j]) - 2 * as.numeric(res.summ[1, j])
    res.summ[4, j] <- as.numeric(res.summ[3, j]) + (2 * as.numeric(res.summ[2, j]) * (as.numeric(res.summ[2, j]) + 1))/(length(div) - as.numeric(res.summ[2, j]) - 1)
    res.summ[7, j] <- res$convergence # termination code from nlminb
    res.summ[8, j] <- res$par[1] # alpha constant
    res.summ[9, j] <- res$par[2] # sig2 constant
    res.summ[10, j] <- res$par[4] # psi1 constant
    res.summ["psi_slope", j] <- res$par[3] # psi1_slope
    res.summ["psi_2", j] <- res$par[6] # psi2 constant
    res.summ["psi2_slope", j] <- res$par[5] # psi2_slope
  }

  # summarize results in output matrix
  #res.summ[res.summ[,"convergence"]>0,]=NA 
  res.summ["AICc",][res.summ["convergence",]>0]=NA # converts AICc to NA for models in which convergence wasn't succesful
  res.summ[5,] = res.summ[4,] - min(res.summ[4,], na.rm=TRUE) # delta AICc calculated from only those models with succesful convergence
  res.summ[6,] = exp(-0.5*res.summ[5,])/sum(exp(-0.5*res.summ[5,]),na.rm=TRUE) # Akaike weights calc'd from deltAICc
  res.summ = apply(res.summ, MARGIN=2, FUN=round, digits=8) 
  return(res.summ)
}
