model_error_rate <-
function(datasets, ages, sim_model, alternatives, type, me1=NULL, me2=NULL, GRAD=NULL, cats=NULL, breakpoint=NULL, domain=NULL, threshold=0, parallel=FALSE, cores=NULL) {
  
  if(sim_model %in% c("BM_null", "OU_null", "BM_linear", "OU_linear", "OU_linear_sig", 
    "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_cat", "DA_OU", 
    "DA_OU_linear", "DA_OU_cat", "DA_BM", "DA_BM_linear", "DA_BM_cat", "OU_BM", "OU_BM_linear", 
    "OU_BM_cat") == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }
  if(is.list(datasets) == FALSE) stop("error rate 'datasets' must be a list with length > 1")
  if(length(datasets) <= 1) stop("error rate 'datasets' must be a list with length > 1")
  if(type != 1 & type != 2) stop("error rate 'type' argument must be set to the integer 1 or 2")

  # perform model selection given the models chosen by the user
  if(parallel == TRUE) {
    if(is.null(cores)) ncor = detectCores()
    if(is.null(cores) == FALSE) ncor = cores
    y = mclapply(X=datasets, FUN=model_select, ages=ages, models=c(sim_model, alternatives), GRAD=GRAD, cats=cats, me1=NULL, me2=NULL, breakpoint=breakpoint, domain=domain, mc.cores=ncor)
  } else {
    y = lapply(X=datasets, FUN=model_select, ages=ages, models=c(sim_model, alternatives), GRAD=GRAD, cats=cats, me1=NULL, me2=NULL, breakpoint=breakpoint, domain=domain)
  }
  
  # calculate error rate
  if(type == 1) {
    #FITS = lapply(y, function(y) y["delta_AICc", sim_model])
    #error_rate = sum(FITS > threshold, na.rm=TRUE)/length(FITS)
    FITS = sapply(y, function(y) sapply(c(sim_model,alternatives), function(k) y["delta_AICc", k]))
    res = apply(FITS, 1, function(x) 1-sum(x > threshold, na.rm=TRUE)/length(x))
    res["err.rate"] = 1 - res[sim_model]
  }
  if(type == 2) {
  	FITS = sapply(y, function(y) sapply(c(sim_model,alternatives), function(k) y["delta_AICc", k]))
  	res = apply(FITS, 1, function(x) 1-sum(x > -threshold, na.rm=TRUE)/length(x))
  	res["err.rate"] = 1 - res[sim_model]
    #FITS = lapply(y, function(y) y["AICc", sim_model] - min(y["AICc", alternatives]))
    #error_rate = sum(FITS > -threshold, na.rm=TRUE)/length(FITS)
  }
  return(res)
}