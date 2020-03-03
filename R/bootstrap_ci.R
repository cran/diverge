bootstrap_ci <-
function(div, ages, GRAD=NULL, cats=NULL, breakpoint=NULL, domain=NULL, model, 
N, parallel=FALSE, cores=NULL, starting=NULL) {

  if(model %in% c("BM_null", "OU_null", "BM_linear", "OU_linear", "OU_linear_sig", 
    "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_cat") == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }

  # Generate random indices
  indices = vector("list", N)
  for (i in 1:N) indices[[i]] = sample(1:length(div), size=length(div), replace=T)

  # re-estimate parameters from bootstrap data
  if(parallel==TRUE) {
    if(is.null(cores)) ncor=detectCores()
    if(is.null(cores)==FALSE) ncor=cores
    boot_pars = mclapply(X=indices, FUN=function(x) find_mle(div=div[x], ages=ages[x], GRAD=GRAD[x], cats=cats[x], domain=domain, breakpoint=breakpoint[x], model=model, p_starting=starting), mc.cores=ncor)
  }
  if(parallel==FALSE) {
    boot_pars = lapply(X=indices, FUN=function(x) find_mle(div=div[x], ages=ages[x], GRAD=GRAD[x], cats=cats[x], domain=domain, breakpoint=breakpoint[x], model=model, p_starting=starting))
  }

  # summarize parameter estimates
  par_matrix = do.call(rbind, lapply(boot_pars, function(x) x$par))
  colnames(par_matrix) = col_names(model=model, ncats=length(unique(cats)))
  par_summ = matrix(nrow = 4, ncol = ncol(par_matrix))
  rownames(par_summ) = c("mean", "median", "percentile_low_95ci", "percentil_high_95ci")
  colnames(par_summ) = colnames(par_matrix)
  for(j in 1:ncol(par_summ)) {
    par_summ["mean", j] = mean(par_matrix[, j])
    par_summ["median", j] = median(par_matrix[, j])
    par_summ[3:4, j] = quantile(x = par_matrix[,j], probs=c(0.025, 0.975))
  }

  # summarize bootstrap info in result list
  res = list(model = model, boostrap_pars = par_matrix, N = N, param_summary = par_summ)
  return(res)
}
