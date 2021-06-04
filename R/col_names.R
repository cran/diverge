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
  if(model == "BM_cat") {
    res=paste("sig2",1:ncats,sep="")
  }
  if(model == "OU_cat") {
    alphas=paste("alpha",1:ncats,sep="")
    res=c("sig2", alphas)
  }
  if(model == "DA_cat") {
    psis=paste("psi",1:ncats,sep="")
    res=c("alpha","sig2",psis)
  }
  if(model == "DA_wt") res=c("alpha", "sig2", "psi1", "psi2", "wt")
  if(model == "DA_bp") res=c("alpha", "sig2", "psi1", "psi2")
  if(model == "DA_wt_linear") res=c("alpha", "sig2", "wt", "psi1_slope", "psi1_int", "psi2_slope", "psi2_int")
  if(model == "DA_bp_linear") res=c("alpha", "sig2", "psi1_slope", "psi1_int", "psi2_slope", "psi2_int")
return(res)
}
