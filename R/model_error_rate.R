model_error_rate <-
function(datasets, ages, sim_model, alternatives, type, GRAD=NULL, cats=NULL, breakpoint=NULL, domain=NULL, threshold=0, parallel=FALSE, cores=NULL) {
  
  if(sim_model %in% c("BM_null", "OU_null", "BM_linear", "OU_linear", "OU_linear_sig", 
    "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_cat") == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }
  if(is.list(datasets) == FALSE) stop("error rate 'datasets' must be a list with length > 1")
  if(length(datasets) <= 1) stop("error rate 'datasets' must be a list with length > 1")
  if(type != 1 & type != 2) stop("error rate 'type' argument must be set to the integer 1 or 2")

  # perform model selection given the models chosen by the user
  if(parallel == TRUE) {
    if(is.null(cores)) ncor = detectCores()
    if(is.null(cores) == FALSE) ncor = cores
    y = mclapply(X=datasets, FUN=model_select, ages=ages, models=c(sim_model, alternatives), GRAD=GRAD, cats=cats, breakpoint=breakpoint, domain=domain, mc.cores=ncor)
  } else {
    y = lapply(X=datasets, FUN=model_select, ages=ages, models=c(sim_model, alternatives), GRAD=GRAD, cats=cats, breakpoint=breakpoint, domain=domain)
  }
  
   # calculate error rate
  if(type == 1) {
    FITS = lapply(y, function(y) y["delta_AICc", sim_model])
    error_rate = sum(FITS > threshold, na.rm=TRUE)/length(FITS)
  }
  if(type == 2) {
    FITS = lapply(y, function(y) y["AICc", sim_model] - min(y["AICc", alternatives]))
    error_rate = sum(FITS > -threshold, na.rm=TRUE)/length(FITS)
  }
  return(error_rate)
}
