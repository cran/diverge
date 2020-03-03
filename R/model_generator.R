model_generator <-
function(domain, par_low, par_high) {
  dom = domain[2]-domain[1]
  slopes = do.call(c, lapply(par_high, function(x) (x-par_low)/dom))
  ints = slopes*(-domain[1]) + par_low
  models=cbind(slopes, ints)
  if(length(par_high)==1) {
    models=models[-which(is.na(models[,1])),]
  }
  if(length(par_low)==1) {
    models=models[-which(is.na(models[,2])),]
  }
  return(models)
}
