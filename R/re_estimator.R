re_estimator <-
function (model, div, ages, me1= NULL, me2 = NULL, GRAD=NULL, cats=NULL, breakpoint=NULL, domain=NULL, p_starting=NULL, parallel=FALSE, cores=NULL, absolute=TRUE) {
  
  if(model %in% c("BM_null", "OU_null", "BM_linear", "OU_linear", "OU_linear_sig", 
    "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_cat") == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }

  if(class(div)=="numeric") { 
    res = find_mle(model=model, p_starting = p_starting, div=div, ages=ages, GRAD=GRAD, cats=cats, breakpoint=breakpoint, domain=domain, absolute=absolute)
    estimate=matrix(res$par, nrow=1)
  }
  if(class(div)=="list") {
    if(parallel==TRUE) {
      if(is.null(cores)) ncor = detectCores()
      if(is.null(cores)==FALSE) ncor = cores
      res = mclapply(div, FUN=find_mle, model = model, p_starting = p_starting, ages = ages, GRAD = GRAD, cats=cats, breakpoint=breakpoint, domain = domain, absolute = absolute, mc.cores=ncor)
      estimate = do.call(rbind, lapply(res, function(x) x$par))
    } else {
      res = lapply(div, FUN=find_mle, model = model, p_starting = p_starting, ages = ages, GRAD = GRAD, cats=cats, breakpoint=breakpoint, domain = domain, absolute = absolute)
      estimate = do.call(rbind, lapply(res, function(x) x$par))
      rownames(estimate)=NULL
    }
  }
  colnames(estimate) = col_names(model, ncats=length(unique(cats)))
  return(estimate)
}

