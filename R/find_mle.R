find_mle <-
function(model, p_starting=NULL, div, ages, GRAD = NULL, cats=NULL, breakpoint=NULL, domain=NULL, absolute=TRUE, parallel=FALSE, cores=NULL) {
    
  if(model %in% c("BM_null", "OU_null", "BM_linear", "OU_linear", "OU_linear_sig", 
    "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_cat") == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }
    # gather a set of starting parameter values
    if (is.null(p_starting)==FALSE) {
      result_matrix <- cbind(p_starting, rep(NA, nrow(p_starting)))
    }
    if (is.null(p_starting)) {
      if(is.null(cats)) {
        result_matrix <- param_grid(model=model, domain=domain)
        } else {
        result_matrix <- param_grid(model=model, domain=domain, ncats=length(unique(cats)))
        }
    }
    param_list = as.list(as.data.frame(t(result_matrix[,1:ncol(result_matrix)-1])))

    # Define boundaries for the parameters in each model
    if (model=="BM_null" | model=="OU_null" | model=="DA_null" | model=="DA_wt" | model=="DA_bp" | model=="DA_cat") lim=1e-5
    if (model=="BM_linear") lim=c(-Inf, 1e-5)
    if (model=="OU_linear") lim=c(1e-5, -Inf, 1e-5)
    if (model=="OU_linear_sig") lim=c(1e-5, -Inf, 1e-5)
    if (model=="DA_linear") lim=c(rep(1e-5, 2), -Inf, 1e-5)
    if (model=="DA_wt_linear") lim=c(rep(1e-5, 3), rep(c(-Inf, 1e-5),2))
    if (model=="DA_bp_linear") lim=c(rep(1e-5,2), rep(c(-Inf, 1e-5),2))
    
    # search likelihood space from the various parameter starting points
    if (parallel==FALSE) {
      res = suppressWarnings(lapply(param_list, FUN=nlminb, objective = likelihood_func, model = model, div = div, ages = ages, cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim, control=list("eval.max"=1000, "iter.max"=5000)))
    }
    if(parallel == TRUE) {
      if(is.null(cores)) {
        ncor=detectCores()
        res <- mclapply(param_list, FUN=nlminb, objective = likelihood_func, model = model, div = div, ages = ages, cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim, control=list("eval.max"=1000, "iter.max"=5000), mc.cores=ncor) 
      } else {
        res = mclapply(param_list, FUN=nlminb, objective = likelihood_func, model = model, div = div, ages = ages, cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim, mc.cores=cores, control=list("eval.max"=1000, "iter.max"=5000))
      }
    }
    # store results in matrix
    result_mat_values = do.call(rbind, lapply(res, function(x) c(x$par, x$objective, x$convergence)))
    
    # set param values and likelihood to NA for the models that didn't converge (i.e. nullify unreliable values)
    result_mat_values[which(result_mat_values[,ncol(result_mat_values)]>0),1:ncol(result_mat_values)-1] = NA 
    
    # set the last column (all NAs) in matrix of starting params to the likelihoods returned from them
    result_matrix[,ncol(result_matrix)] = result_mat_values[,ncol(result_mat_values)-1] 
    
    # order the matrices in increasing order by neg.log.likelihood
    result_mat_values_ordered = result_mat_values[order(as.numeric(result_mat_values[,ncol(result_mat_values)-1])),] 
    result_mat_ordered = result_matrix[order(as.numeric(result_matrix[,ncol(result_matrix)])),]
    
    # re-run the likelihood search from a starting parameter value from which the likelihood search returned the mle
    # if the search doesn't return the same end values, it was stuck on a local optima and the final estimate is not the true mle
    p = as.numeric(result_mat_ordered[1, 1:(ncol(result_mat_ordered)-1)])
    res = suppressWarnings(nlminb(start=p, objective = likelihood_func, model = model, div = div, ages = ages, cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim, control=list("eval.max"=1000, "iter.max"=5000)))
    if (round(res$objective, 4) > round(result_mat_values_ordered[1, ncol(result_mat_values_ordered)-1], 4)) {
       warning(paste(model,": final result is not MLE", sep=""))
    }
    names(res$par) = col_names(model, ncats=length(unique(cats)))
  return(res)
}

