param_grid <-
function (model, domain=NULL, ncats=NULL) {

  if(model %in% c("BM_null", "OU_null", "BM_linear", "OU_linear", "OU_linear_sig", 
    "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_cat") == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }
   
    if (model == "BM_null") {
      beta = c(0.001, 0.01, 0.1, 1, 10, 100)
      result_matrix = cbind(beta, rep(NA, length(beta)))
    }
    if (model == "BM_linear") {
      if(is.null(domain)==TRUE) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      beta_start = c(0.001, 0.01, 0.1, 1, 10, 20, 50, 100)
      beta_models = model_generator(par_low=beta_start, par_high=beta_start, domain=domain)
      result_matrix = cbind(beta_models, rep(NA, nrow(beta_models)))
    }
    if (model == "OU_null") {
      alpha = c(0.001, 0.01, 0.1, 1, 10, 100)
      beta = c(0.001, 0.01, 0.1, 1, 10, 100)
      paramgrid=as.matrix(expand.grid(alpha,beta))
      result_matrix=cbind(paramgrid, rep(NA, nrow(paramgrid)))
    }
    if(model == "OU_linear") {
      if(is.null(domain)) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      beta = c(0.001, 0.01, 0.1, 1, 10, 100)
      alpha_start = c(0.001, 0.01, 0.1, 1, 10, 100)
      alpha_models = model_generator(par_low=alpha_start, par_high=alpha_start, domain=domain)
      result_matrix = matrix(NA, nrow=nrow(alpha_models)*length(beta), ncol = 4)
      counter=0
      for(i in 1:nrow(alpha_models)) {
        models=alpha_models[i,]
          for(k in 1:length(beta)) {
            counter=counter+1
            result_matrix[counter, 1] = models[2]
            result_matrix[counter, 2] = models[1]
            result_matrix[counter, 3] = beta[k]
          }
      }
    }
    if(model == "OU_linear_sig") {
      if(is.null(domain)) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      alpha = c(0.001, 0.01, 0.1, 1, 10, 100)
      beta_start = c(0.001, 0.01, 0.1, 1, 10, 100)
      beta_models = model_generator(par_low=beta_start, par_high=beta_start, domain=domain)
      result_matrix = matrix(NA, nrow=nrow(beta_models)*length(alpha), ncol = 4)
      counter=0
      for(i in 1:nrow(beta_models)) {
        models=beta_models[i,]
          for(k in 1:length(alpha)) {
            counter=counter+1
            result_matrix[counter, 1] = alpha[k]
            result_matrix[counter, 2] = models[1]
            result_matrix[counter, 3] = models[2]
          }
      }
    }
    if (model == "DA_null") {
      alpha = c(0.001, 0.01, 0.1, 1, 10, 100)
      beta = c(0.001, 0.01, 0.1, 1, 10, 100)
      psi = seq(0.0001,2,0.2)
      paramgrid=as.matrix(expand.grid(alpha, beta, psi))
      result_matrix=cbind(paramgrid, rep(NA, nrow(paramgrid)))
    }
    if (model == "DA_linear") {
      if(is.null(domain)) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      alpha = c(0.001, 0.01, 0.1, 1, 10, 100)
      beta = c(0.001, 0.01, 0.1, 1, 10, 100)
      rates = as.matrix(expand.grid(alpha, beta))
      psi_start= c(0.0005, 0.05, 0.5, 5, 50)
      psi_end=psi_start
      psi_models = model_generator(par_low=psi_start, par_high=psi_end, domain=domain) 
      result_matrix = matrix(NA, nrow=nrow(rates)*nrow(psi_models), ncol = 5)
      counter=0
      for(i in 1:nrow(rates)) {
          rate_pars=rates[i,]
          for(k in 1:nrow(psi_models)) {
              models=psi_models[k,]
              counter=counter+1
              result_matrix[counter, 1:2] <- rate_pars
              result_matrix[counter, 3:4] <- models
            }
        }
    }
    if (model == "DA_cat") {
      if(is.null(ncats)) stop("Number of categories in D_cat not detected")
      if(ncats > 3) stop("More than the maximum of 3 categories have been provided")
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      beta = c(0.001, 0.01, 0.1, 1, 10)
      psi1 = c(0.0005, 0.05, 0.5, 5, 50)
      psi2 = psi1
      if(ncats == 2) paramgrid=as.matrix(expand.grid(alpha, beta, psi1, psi2))
      if(ncats == 3) {
        psi3 = psi1
        paramgrid=as.matrix(expand.grid(alpha, beta, psi1, psi2, psi3))
      }
      result_matrix=cbind(paramgrid, rep(NA, length=nrow(paramgrid)))
    }
    if (model == "DA_wt") {
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      beta = c(0.001, 0.01, 0.1, 1, 10)
      psi1 = seq(0.0001,2,0.2)
      psi2 = psi1
      bp = c(0.001, seq(2, 10, 2))
      paramgrid=as.matrix(expand.grid(alpha, beta, psi1, psi2, bp))
      result_matrix=cbind(paramgrid, matrix(NA, nrow=nrow(paramgrid), ncol=1))
    }
    if (model == "DA_bp") {
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      beta = c(0.001, 0.01, 0.1, 1, 10)
      psi1 = seq(0.0001,2,0.2)
      psi2 = psi1
      paramgrid=as.matrix(expand.grid(alpha, beta, psi1, psi2))
      result_matrix=cbind(paramgrid, matrix(NA, nrow=nrow(paramgrid), ncol=1))
    }
    if (model == "DA_wt_linear") {
      if(is.null(domain)) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      beta = c(0.001, 0.01, 0.1, 1)
      wt = c(0.001, 1, 10)
      rates = as.matrix(expand.grid(alpha, beta, wt))
      psi_start= c(0.0005, 0.05, 0.5, 5, 50)
      psi_end=psi_start
      psi_models = model_generator(par_low=psi_start, par_high=psi_end, domain=domain)
      psi2_models = psi_models
      result_matrix = matrix(NA, nrow=nrow(rates)*nrow(psi_models)*nrow(psi_models), ncol = 8)
      counter=0
      for(i in 1:nrow(rates)) {
          rate_pars=rates[i,]
          for(k in 1:nrow(psi_models)) {
              models=psi_models[k,]
              for(j in 1:nrow(psi2_models)) {
                models2 = psi2_models[j,]
                counter=counter+1
                result_matrix[counter, 1:3] <- rate_pars
                result_matrix[counter, 4:5] <- models
                result_matrix[counter, 6:7] <- models2
            }
          }
        }
    }
    if (model == "DA_bp_linear") {
      if(is.null(domain)) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      beta = c(0.001, 0.01, 0.1, 1)
      rates = as.matrix(expand.grid(alpha, beta))
      psi_start= c(0.0005, 0.05, 0.5, 5, 50)
      psi_end=psi_start
      psi_models = model_generator(par_low=psi_start, par_high=psi_end, domain=domain)
      psi2_models = psi_models
      result_matrix = matrix(NA, nrow=nrow(rates)*nrow(psi_models)*nrow(psi_models), ncol = 7)
      counter=0
      for(i in 1:nrow(rates)) {
          rate_pars=rates[i,]
          for(k in 1:nrow(psi_models)) {
              models=psi_models[k,]
              for(j in 1:nrow(psi2_models)) {
                models2 = psi2_models[j,]
                counter=counter+1
                result_matrix[counter, 1:2] <- rate_pars
                result_matrix[counter, 3:4] <- models
                result_matrix[counter, 5:6] <- models2
            }
          }
        }
    }
 colnames(result_matrix) = c(col_names(model, ncats=ncats),"likelihood")
 return(result_matrix)
}