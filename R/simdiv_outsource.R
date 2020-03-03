simdiv_outsource <-
function(u2, var, model) {
  if(model != "BM_null" & model != "OU_null") {
    if(length(u2) != length(var)) {
      stop("your expected value and variance vectors are of different lengths")
    }
  }  
	sis1=rnorm(n=length(var), mean=u2, sd=sqrt(var))
	sis2=rnorm(n=length(var), mean=0, sd=sqrt(var))
	div_sis=abs(sis2-sis1)
  return(div_sis)
}
