param_grid <-
function (model, domain=NULL, ncats=NULL) {

  if(model %in% c("BM_null", "BM_linear", "BM_cat", "OU_null", "OU_linear", "OU_linear_sig", "OU_cat",
    "DA_null", "DA_linear", "DA_cat", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_OU", "DA_BM",
    "DA_OU_linear", "DA_BM_linear", "DA_OU_cat", "DA_BM_cat", "OU_BM", "OU_BM_linear", "OU_BM_cat") == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }
   
    if (model == "BM_null") {
      sig2 = c(0.001, 0.01, 0.1, 1, 10, 100)
      result_matrix = cbind(sig2, "negloglike"=rep(NA, length(sig2)))
    }
    if (model == "BM_linear") {
      if(is.null(domain)==TRUE) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      sig2_start = c(0.001, 0.01, 0.1, 1, 10, 20, 50, 100)
      sig2_models = model_generator(par_low=sig2_start, par_high=sig2_start, domain=domain)
      colnames(sig2_models) = c("sig2_slope", "sig2_int")
      result_matrix = cbind(sig2_models, "negloglike"=rep(NA, nrow(sig2_models)))
    }
    if (model == "OU_null") {
      alpha = c(0.001, 0.01, 0.1, 1, 10, 100)
      sig2 = c(0.001, 0.01, 0.1, 1, 10, 100)
      paramgrid=as.matrix(expand.grid(alpha=alpha,sig2=sig2))
      result_matrix=cbind(paramgrid, "negloglike"=rep(NA, nrow(paramgrid)))
    }
    if(model == "OU_linear") {
      if(is.null(domain)) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      sig2 = c(0.001, 0.01, 0.1, 1, 10, 100)
      alpha_start = c(0.001, 0.01, 0.1, 1, 10, 100)
      alpha_models = model_generator(par_low=alpha_start, par_high=alpha_start, domain=domain)
      result_matrix = matrix(NA, nrow=nrow(alpha_models)*length(sig2), ncol = 4)
      colnames(result_matrix) = c("alpha_int", "alpha_slope", "sig2", "negloglike")
      counter=0
      for(i in 1:nrow(alpha_models)) {
        models=alpha_models[i,]
          for(k in 1:length(sig2)) {
            counter=counter+1
            result_matrix[counter, 1] = models[2]
            result_matrix[counter, 2] = models[1]
            result_matrix[counter, 3] = sig2[k]
          }
      }
    }
    if(model == "OU_linear_sig") {
      if(is.null(domain)) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      alpha = c(0.001, 0.01, 0.1, 1, 10, 100)
      sig2_start = c(0.001, 0.01, 0.1, 1, 10, 100)
      sig2_models = model_generator(par_low=sig2_start, par_high=sig2_start, domain=domain)
      result_matrix = matrix(NA, nrow=nrow(sig2_models)*length(alpha), ncol = 4)
      colnames(result_matrix) = c("alpha_int", "sig2_slope", "sig2_int", "negloglike")
      counter=0
      for(i in 1:nrow(sig2_models)) {
        models=sig2_models[i,]
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
      sig2 = c(0.001, 0.01, 0.1, 1, 10, 100)
      psi = seq(0.0001,2,0.2)
      paramgrid=as.matrix(expand.grid(alpha=alpha, sig2=sig2, psi=psi))
      result_matrix=cbind(paramgrid, "negloglike"=rep(NA, nrow(paramgrid)))
    }
    if (model == "DA_linear") {
      if(is.null(domain)) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      alpha = c(0.001, 0.01, 0.1, 1, 10, 100)
      sig2 = c(0.001, 0.01, 0.1, 1, 10, 100)
      rates = as.matrix(expand.grid(alpha, sig2))
      psi_start= c(0.0005, 0.05, 0.5, 5, 50)
      psi_end=psi_start
      psi_models = model_generator(par_low=psi_start, par_high=psi_end, domain=domain) 
      result_matrix = matrix(NA, nrow=nrow(rates)*nrow(psi_models), ncol = 5)
      colnames(result_matrix) = c("alpha", "sig2", "psi_slope", "psi_int", "negloglike")
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
      if(is.null(ncats)) stop("Number of categories not detected")
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      psi = c(0.0005, 0.05, 0.5, 5, 50)
      psidf = as.data.frame(sapply(1:ncats, function(x) psi))
      for(i in 1:ncats) colnames(psidf)[i] = paste("psi",i,sep="")
      pars = cbind(alpha, sig2, psidf)
      paramgrid = as.matrix(expand.grid(pars))
      result_matrix=cbind(paramgrid, "negloglike"=rep(NA, length=nrow(paramgrid)))
    }
    if (model == "BM_cat") {
      if(is.null(ncats)) stop("Number of categories not detected")
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      sigdf = as.data.frame(sapply(1:ncats, function(x) sig2))
      for(i in 1:ncats) colnames(sigdf)[i] = paste("sig2_",i,sep="")
      paramgrid = as.matrix(expand.grid(sigdf))
      result_matrix=cbind(paramgrid, "negloglike"=rep(NA, length=nrow(paramgrid)))
    }
    if (model == "OU_cat") {
      if(is.null(ncats)) stop("Number of categories not detected")
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      alphadf = as.data.frame(sapply(1:ncats, function(x) alpha))
      for(i in 1:ncats) colnames(alphadf)[i] = paste("alpha",i,sep="")
      pars=cbind(sig2,alphadf)
      paramgrid = as.matrix(expand.grid(pars))
      result_matrix=cbind(paramgrid, "negloglike"=rep(NA, length=nrow(paramgrid)))
    }
    if (model == "DA_OU" | model == "DA_BM") {
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      psi = c(0.0005, 0.05, 0.5, 5, 50)
      prop = c(0.1, 0.5, 0.9)
      paramgrid = as.matrix(expand.grid(alpha=alpha, sig2=sig2, psi=psi, prop=prop))
      result_matrix=cbind(paramgrid, "negloglike"=rep(NA, length=nrow(paramgrid)))
    }
    if (model == "OU_BM") {
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      prop = c(0.1, 0.5, 0.9)
      paramgrid = as.matrix(expand.grid(alpha=alpha, sig2=sig2, prop=prop))
      result_matrix=cbind(paramgrid, "negloglike"=rep(NA, length=nrow(paramgrid)))
    }
    if (model == "DA_OU_linear" | model == "DA_BM_linear") {
      if(is.null(domain)) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      psi = c(0.005, 0.05, 0.5, 5, 50)
      pars = expand.grid(alpha, sig2, psi)
      prop_start = prop_end = c(0.1, 0.5, 0.9)
      prop_models = model_generator(par_low = prop_start, par_high=prop_end, domain=domain)
      result_matrix = matrix(NA, nrow=nrow(pars)*nrow(prop_models), ncol = 6)
      colnames(result_matrix) = c("alpha", "sig2", "psi", "prop_slope", "prop_int", "negloglike")
      counter=0
      for(i in 1:nrow(pars)) {
          rate_pars=as.numeric(pars[i,])
          for(k in 1:nrow(prop_models)) {
              models=as.numeric(prop_models[k,])
              counter=counter+1
              result_matrix[counter, 1:3] <- rate_pars
              result_matrix[counter, 4:5] <- models
          }
      }
    }
    if (model == "OU_BM_linear") {
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      pars = expand.grid(alpha, sig2)
      prop_start = prop_end = c(0.1, 0.5, 0.9)
      prop_models = model_generator(par_low = prop_start, par_high=prop_end, domain=domain)
      result_matrix = matrix(NA, nrow=nrow(pars)*nrow(prop_models), ncol = 5)
      colnames(result_matrix) = c("alpha", "sig2", "prop_slope", "prop_int", "negloglike")
      counter=0
      for(i in 1:nrow(pars)) {
          rate_pars=as.numeric(pars[i,])
          for(k in 1:nrow(prop_models)) {
              models=as.numeric(prop_models[k,])
              counter=counter+1
              result_matrix[counter, 1:2] <- rate_pars
              result_matrix[counter, 3:4] <- models
          }
      }
    }
    if (model == "DA_OU_cat" | model == "DA_BM_cat") {
      if(is.null(ncats)) stop("Number of categories not detected")
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      psi = c(0.005, 0.05, 0.5, 5, 50)
      pars=as.matrix(expand.grid(alpha, sig2, psi))
      prop = c(0.1, 0.5, 0.9)
      propdf = as.data.frame(sapply(1:ncats, function(x) prop))
      propdf = as.matrix(expand.grid(propdf))
      for(i in 1:ncats) colnames(propdf)[i] = paste("prop",i,sep="")
      result_matrix = matrix(NA, nrow=nrow(pars)*nrow(propdf), ncol = 4+ncats)
      colnames(result_matrix) = c("alpha", "sig2", "psi", paste("prop", 1:ncats, sep=""),"negloglike")
      counter=0
      for(i in 1:nrow(pars)) {
          rate_pars=as.numeric(pars[i,])
          for(k in 1:nrow(propdf)) {
              propns=as.numeric(propdf[k,])
              x = 3+ncats
              counter=counter+1
              result_matrix[counter, 1:3] <- rate_pars
              result_matrix[counter, 4:x] <- propns
          }
      }
    }
    if (model == "OU_BM_cat") {
      if(is.null(ncats)) stop("Number of categories not detected")
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      pars=as.matrix(expand.grid(alpha, sig2))
      prop = c(0.1, 0.5, 0.9)
      propdf = as.data.frame(sapply(1:ncats, function(x) prop))
      propdf = as.matrix(expand.grid(propdf))
      for(i in 1:ncats) colnames(propdf)[i] = paste("prop",i,sep="")
      result_matrix = matrix(NA, nrow=nrow(pars)*nrow(propdf), ncol = 3+ncats)
      colnames(result_matrix) = c("alpha", "sig2", paste("prop", 1:ncats, sep=""),"negloglike")
      counter=0
      for(i in 1:nrow(pars)) {
          rate_pars=as.numeric(pars[i,])
          for(k in 1:nrow(propdf)) {
              propns=as.numeric(propdf[k,])
              x = 2+ncats
              counter=counter+1
              result_matrix[counter, 1:2] <- rate_pars
              result_matrix[counter, 3:x] <- propns
          }
      }
    }
    if (model == "DA_wt") {
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      psi1 = seq(0.0001,2,0.2)
      psi2 = psi1
      wt = c(0.001, seq(2, 10, 2))
      paramgrid=as.matrix(expand.grid(alpha=alpha, sig2=sig2, psi1=psi1, psi2=psi2, wt=wt))
      result_matrix=cbind(paramgrid, "negloglike"=rep(NA, length=nrow(paramgrid)))
    }
    if (model == "DA_bp") {
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1, 10)
      psi1 = seq(0.0001,2,0.2)
      psi2 = psi1
      paramgrid=as.matrix(expand.grid(alpha=alpha, sig2=sig2, psi1=psi1, psi2=psi2))
      result_matrix=cbind(paramgrid, "negloglike"=rep(NA, length=nrow(paramgrid)))
    }
    if (model == "DA_wt_linear") {
      if(is.null(domain)) {
        stop("Hold up: You've chosen a model that requires a domain argument")
      }
      alpha = c(0.001, 0.01, 0.1, 1, 10)
      sig2 = c(0.001, 0.01, 0.1, 1)
      wt = c(0.001, 1, 10)
      rates = as.matrix(expand.grid(alpha, sig2, wt))
      psi_start= c(0.0005, 0.05, 0.5, 5, 50)
      psi_end=psi_start
      psi_models = model_generator(par_low=psi_start, par_high=psi_end, domain=domain)
      psi2_models = psi_models
      result_matrix = matrix(NA, nrow=nrow(rates)*nrow(psi_models)*nrow(psi_models), ncol = 8)
      colnames(result_matrix) = c("alpha", "sig2", "wt", "psi1_slope", "psi1_int", 
        "psi2_slope", "psi2_int", "negloglike")
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
      sig2 = c(0.001, 0.01, 0.1, 1)
      rates = as.matrix(expand.grid(alpha, sig2))
      psi_start= c(0.0005, 0.05, 0.5, 5, 50)
      psi_end=psi_start
      psi_models = model_generator(par_low=psi_start, par_high=psi_end, domain=domain)
      psi2_models = psi_models
      result_matrix = matrix(NA, nrow=nrow(rates)*nrow(psi_models)*nrow(psi_models), ncol = 7)
      colnames(result_matrix) = c("alpha", "sig2", "psi1_slope", "psi1_int", 
        "psi2_slope", "psi2_int", "negloglike")
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
 return(result_matrix)
}
