find_mle <-
function(model, p_starting=NULL, div, ages, me1 = NULL, me2 = NULL, GRAD = NULL, cats=NULL, 
    breakpoint=NULL, domain=NULL, absolute=TRUE, parallel=FALSE, cores=NULL) {
    
  if(model %in% c("BM_null", "BM_linear", "BM_cat", "OU_null", "OU_linear", "OU_linear_sig", "OU_cat",
    "DA_null", "DA_linear", "DA_cat", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_OU", "DA_OU_linear",
    "DA_OU_cat", "DA_BM", "DA_BM_linear", "DA_BM_cat", "OU_BM", "OU_BM_linear", "OU_BM_cat") == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }
  if(sum(is.null(me1), is.null(me2)) == 1) stop("You've supplied measurement error 
    for just one species in the pairs; please provide 2 or 0 vectors for measurement error")
    
    # gather a set of starting parameter values
    if (is.null(p_starting)==FALSE) {
      par_matrix <- cbind(p_starting, rep(NA, nrow(p_starting)))
    }
    if (is.null(p_starting)) {
      if(is.null(cats)) {
        par_matrix <- param_grid(model=model, domain=domain)
        } else {
        par_matrix <- param_grid(model=model, ncats=length(unique(cats)))
        }
    }
    param_list = as.list(as.data.frame(t(par_matrix[,1:ncol(par_matrix)-1])))

    # Define boundaries for the parameters in each model
    if (model=="BM_null" | model=="OU_null" | model=="DA_null" | model=="DA_wt" | model=="DA_bp" | 
        model=="DA_cat" | model=="BM_cat" | model=="OU_cat") lim=1e-8
    if (model == "DA_OU" | model == "DA_BM") lim=c(rep(1e-8,3),2/length(div), c(rep(Inf,3),1))
    if (model == "OU_BM") lim=c(rep(1e-8, 2), 2/length(div), rep(Inf, 2), 1)
    if (model == "DA_OU_linear" | model == "DA_BM_linear") lim=c(rep(1e-8,3),-Inf, 2/length(div), c(rep(Inf,4),1))
    if (model == "OU_BM_linear") lim=c(rep(1e-8,2),-Inf, 2/length(div), c(rep(Inf,3),1))
    if (model == "DA_OU_cat" | model == "DA_BM_cat") lim=c(rep(1e-8,3), rep(0,length(unique(cats))), rep(Inf, 3), rep(1, length(unique(cats))))
    if (model == "OU_BM_cat") lim=c(rep(1e-8,2), rep(0,length(unique(cats))), rep(Inf, 2), rep(1, length(unique(cats))))
    if (model=="BM_linear") lim=c(-Inf, 1e-8)
    if (model=="OU_linear") lim=c(1e-8, -Inf, 1e-8)
    if (model=="OU_linear_sig") lim=c(1e-8, -Inf, 1e-8)
    if (model=="DA_linear") lim=c(rep(1e-8, 2), -Inf, 1e-8)
    if (model=="DA_wt_linear") lim=c(rep(1e-8, 3), rep(c(-Inf, 1e-8),2))
    if (model=="DA_bp_linear") lim=c(rep(1e-8,2), rep(c(-Inf, 1e-8),2))
    
    # search likelihood space from the various parameter starting points
    if (parallel==FALSE) {
      if(model=="DA_OU" | model=="DA_BM" | model=="OU_BM" | model=="DA_OU_linear" | model=="DA_BM_linear" | model=="OU_BM_linear" |
        model=="DA_OU_cat" | model=="DA_BM_cat" | model=="OU_BM_cat") {
        res = suppressWarnings(lapply(param_list, FUN=nlminb, objective = likelihood_func, model = model, div = div, ages = ages, me1 = me1, me2 = me2, cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim[1:(length(lim)/2)], upper=lim[(length(lim)/2+1):length(lim)], control=list("eval.max"=1000, "iter.max"=5000)))
      } else {
        res = suppressWarnings(lapply(param_list, FUN=nlminb, objective = likelihood_func, model = model, div = div, ages = ages, me1 = me1, me2 = me2, cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim, control=list("eval.max"=1000, "iter.max"=5000)))
      #res = lapply(param_list, FUN=nlminb, objective = likelihood_func, model = model, div = div, ages = ages, me1 = me1, me2 = me2, cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim, control=list("eval.max"=1000, "iter.max"=5000))
      }
    }
    if(parallel == TRUE) {
      if(is.null(cores)) {
        ncor=detectCores()
        if(model=="DA_OU" | model=="DA_BM" | model=="OU_BM" | model=="DA_OU_linear" | model=="DA_BM_linear" | model=="OU_BM_linear" |
         model=="DA_OU_cat" | model=="DA_BM_cat" | model=="OU_BM_cat") {
          res <- mclapply(param_list, FUN=nlminb, objective = likelihood_func, model = model, div = div, ages = ages, me1 = me1, me2 = me2,cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim[1:(length(lim)/2)], upper=lim[(length(lim)/2+1):length(lim)], control=list("eval.max"=1000, "iter.max"=5000), mc.cores=ncor) 
        } else {
          res <- mclapply(param_list, FUN=nlminb, objective = likelihood_func, model = model, div = div, ages = ages, me1 = me1, me2 = me2,cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim, control=list("eval.max"=1000, "iter.max"=5000), mc.cores=ncor) 
        }
      } else {
        if(model=="DA_OU" | model=="DA_BM" | model=="OU_BM" | model=="DA_OU_linear" | model=="DA_BM_linear" | model=="OU_BM_linear" |
         model=="DA_OU_cat" | model=="DA_BM_cat" | model=="OU_BM_cat") {
          res = mclapply(param_list, FUN=nlminb, objective = likelihood_func, model = model, div = div, ages = ages, me1 = me1, me2 = me2,cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim[1:(length(lim)/2)], upper=lim[(length(lim)/2+1):length(lim)], mc.cores=cores, control=list("eval.max"=1000, "iter.max"=5000)) 
        } else {
          res = mclapply(param_list, FUN=nlminb, objective = likelihood_func, model = model, div = div, ages = ages, me1 = me1, me2 = me2,cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim, mc.cores=cores, control=list("eval.max"=1000, "iter.max"=5000))
        }
      }
    }

    # store results in a matrix
    results = do.call(rbind, lapply(res, function(x) c(x$par, x$objective, x$convergence)))
    
    # set param values and likelihood to NA for the models that didn't converge (i.e. nullify unreliable values)
    # potential problem here -- sometimes failed fits are returning a 'zero' for convergence and a 'favourable' return message -- an issue with nlminb?
    results[which(results[,ncol(results)]>0),1:(ncol(results)-1)] = NA
    
    # set the last column (all NAs) in matrix of starting params to the likelihoods returned from them
    par_matrix[,ncol(par_matrix)] = results[,ncol(results)-1] 
    
    # order the matrices in increasing order by neg.log.likelihood
    # note that nlminb is a minimizing algorithm, and it is minimizing a function for the negative log likelihood (coded in likelihood_func)
    # thus, smaller values correspond to larger likelihoods and better model fits
    results_ordered = par_matrix[order(as.numeric(par_matrix[,ncol(par_matrix)])),]
    
    # re-run the likelihood search from a starting parameter value from which the likelihood search returned the mle
    # if the search doesn't return the same end values, it was stuck on a local optima and the final estimate is not the true mle
    p = as.numeric(results_ordered[1, 1:(ncol(results_ordered)-1)])
    if(model=="DA_OU" | model=="DA_BM" | model=="OU_BM" | model=="DA_OU_linear" | model=="DA_BM_linear" | model=="OU_BM_linear" |
      model=="DA_OU_cat" | model=="DA_BM_cat" | model=="OU_BM_cat") {
      res = suppressWarnings(nlminb(start=p, objective = likelihood_func, model = model, div = div, ages = ages, me1 = me1, me2 = me2, cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim[1:(length(lim)/2)], upper=lim[(length(lim)/2+1):length(lim)], control=list("eval.max"=1000, "iter.max"=5000)))
    } else {
      res = suppressWarnings(nlminb(start=p, objective = likelihood_func, model = model, div = div, ages = ages, me1 = me1, me2 = me2, cats=cats, GRAD=GRAD, bp=breakpoint, absolute=absolute, lower=lim, control=list("eval.max"=1000, "iter.max"=5000)))
    }
    if (all(results[,ncol(results)]>0) | round(res$objective, 4) > round(results_ordered[1, ncol(results_ordered)], 4)) {
       warning(paste(model,": final result is not MLE", sep=""))
    }
    names(res$par) = colnames(par_matrix)[-length(colnames(par_matrix))]
  return(res)
}