col_names <-
function(model, ncats=NULL) {
  if(model=="DA_cat" & is.null(ncats)) stop("No. categories must be provided in 'cats' argument")
  if(model == "BM_null") res="sig2"
  if(model == "OU_null") res=c("alpha","sig2")
  if(model == "DA_null") res=c("alpha", "sig2", "psi")
  if(model == "BM_linear") res=c("sig2_slope", "sig2_int")
  if(model == "OU_linear") res=c("alpha_int", "alpha_slope", "sig2")
  if(model == "OU_linear_sig") res=c("alpha", "sig2_slope", "sig2_int")
  if(model == "DA_linear") res=c("alpha", "sig2", "psi_slope", "psi_int")
  if(is.null(ncats)==FALSE) {
    if(model == "DA_cat" & ncats==2) res=c("alpha", "sig2", "psi1", "psi2")
    if(model == "DA_cat" & ncats==3) res=c("alpha", "sig2", "psi1", "psi2", "psi3")
  }
  if(model == "DA_wt") res=c("alpha", "sig2", "psi1", "psi2", "wt")
  if(model == "DA_bp") res=c("alpha", "sig2", "psi1", "psi2")
  if(model == "DA_wt_linear") res=c("alpha", "sig2", "wt", "psi1_slope", "psi1_int", "psi2_slope", "psi2_int")
  if(model == "DA_bp_linear") res=c("alpha", "sig2", "psi1_slope", "psi1_int", "psi2_slope", "psi2_int")
return(res)
}
