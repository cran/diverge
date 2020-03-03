random_walks = 
function (model, TIME, nsim, sig2, alpha = NULL, psi = NULL, 
    psi2 = NULL, wt = NULL, theta = 0, centre = 0, steps=100, plot = TRUE, 
    col=c("black","red"), labels = TRUE, ylim = NULL, ...) 
{

  if(round(nsim/2) != nsim/2) stop("Please provide an even number for 'nsim'")
    if (model != "BM_null" & model != "OU_null" & model != "DA_null" & 
        model != "DA_wt" & model != "character_displacement" & 
        model != "sorting" & model != "time_only") {
        stop("Spell check: the model you've entered doesn't match the models accepted by this function")
    }
    if (model == "OU_null" & is.null(alpha)) {
        stop("You've chosen an OU model and haven't provided a value for alpha")
    }
    if (model == "DA_null" & is.null(psi)) {
        stop("You've chosen a 'DA' model and haven't provided a value for psi")
    }
    if (model == "DA_wt" & is.null(wt)) {
        stop("You've chosen the 'DA-wait-time' model and haven't provided a value for wt")
    }
    if (model == "character_displacement" & is.null(wt)) {
      stop("You've chosen the character displacement model and haven't provided a value for wt (i.e. time-to-sympatry)")
    }
    steps = TIME*steps
    dt = TIME/steps
    if (model == "BM_null") {
        dw = matrix(rnorm(nsim*steps, 0, sqrt(dt)), nrow = nsim, 
            ncol = steps)
        y = matrix(0, nrow = nsim, ncol = steps)
        for (k in 1:nrow(y)) {
            for (i in 2:ncol(y)) {
                y[k, i] <- y[k, i - 1] + sqrt(sig2) * dw[k, i - 1]
            }
        }
    }
    if (model == "OU_null") {
        dw = matrix(rnorm(nsim*steps, 0, sqrt(dt)), nrow = nsim, 
            ncol = steps)
        y = matrix(0, nrow = nsim, ncol = steps)
        for (k in 1:nrow(y)) {
            for (i in 2:ncol(y)) {
                y[k, i] <- y[k, i - 1] + alpha * (theta - y[k, 
                  i - 1]) * dt + sqrt(sig2) * dw[k, i - 1]
            }
        }
    }
    if (model == "DA_null") {
        theta1 = centre + (psi/2)
        theta2 = centre - (psi/2)
        dw = matrix(rnorm(nsim * steps/2, 0, sqrt(dt)), nrow = nsim/2, 
            ncol = steps)
        y = matrix(0, nrow = nsim/2, ncol = steps)
        for (k in 1:nrow(y)) {
            for (i in 2:ncol(y)) {
                y[k, i] <- y[k, i - 1] + alpha * (theta1 - y[k, 
                  i - 1]) * dt + sqrt(sig2) * dw[k, i - 1]
            }
        }
        dw2 = matrix(rnorm(nsim * steps/2, 0, sqrt(dt)), nrow = nsim/2, 
            ncol = steps)
        y2 = matrix(0, nrow = nsim/2, ncol = steps)
        for (j in 1:nrow(y2)) {
            for (l in 2:ncol(y2)) {
                y2[j, l] <- y2[j, l - 1] + alpha * (theta2 - 
                  y2[j, l - 1]) * dt + sqrt(sig2) * dw2[j, l - 1]
            }
        }
        y = list(y, y2)
    }
    if (model == "DA_wt") {
      bk = wt*steps/TIME
        theta1 = centre - (psi/2)
        theta2 = centre - (psi2/2)
        theta3 = centre + (psi/2)
        theta4 = centre + (psi2/2)
        dw = matrix(rnorm(nsim * steps/2, 0, sqrt(dt)), nrow = nsim/2, 
            ncol = steps)
        y = matrix(0, nrow = nsim/2, ncol = steps)
        for (k in 1:nrow(y)) {
            for (i in 2:bk) {
                y[k, i] <- y[k, i - 1] + alpha * (theta3 - y[k, 
                  i - 1]) * dt + sqrt(sig2) * dw[k, i - 1]
            }
            for (j in (bk + 1):ncol(y)) {
                y[k, j] <- y[k, j - 1] + alpha * (theta4 - y[k, 
                  j - 1]) * dt + sqrt(sig2) * dw[k, j - 1]
            }
        }
        dw2 = matrix(rnorm(nsim * steps/2, 0, sqrt(dt)), nrow = nsim/2, 
            ncol = steps)
        y2 = matrix(0, nrow = nsim/2, ncol = steps)
        for (j in 1:nrow(y2)) {
            for (l in 2:bk) {
                y2[j, l] <- y2[j, l - 1] + alpha * (theta1 - 
                  y2[j, l - 1]) * dt + sqrt(sig2) * dw2[j, l - 1]
            }
            for (m in (bk + 1):ncol(y2)) {
                y2[j, m] <- y2[j, m - 1] + alpha * (theta2 - 
                  y2[j, m - 1]) * dt + sqrt(sig2) * dw2[j, m - 1]
            }
        }
        y = list(y, y2)
    }
    if (model == "character_displacement") {
      if(wt > TIME) stop("Hold up! Wait time (wt) must be less than total time (TIME)")
      bk = wt*steps/TIME
        theta1 = centre - (psi/2)
        theta2 = centre - (psi2/2)
        theta3 = centre + (psi/2)
        theta4 = centre + (psi2/2)
        # matrices for sympatric pairs
        dw_sym1 = matrix(rnorm(nsim * steps/2, 0, sqrt(dt)), nrow = nsim/2, 
            ncol = steps)
        y_sym1 = matrix(0, nrow = nsim/2, ncol = steps)
        for (k in 1:nrow(y_sym1)) {
            for (i in 2:bk) {
                y_sym1[k, i] <- y_sym1[k, i - 1] + alpha * (theta3 - y_sym1[k, 
                  i - 1]) * dt + sqrt(sig2) * dw_sym1[k, i - 1]
            }
            for (j in (bk + 1):ncol(y_sym1)) {
                y_sym1[k, j] <- y_sym1[k, j - 1] + alpha * (theta4 - y_sym1[k, 
                  j - 1]) * dt + sqrt(sig2) * dw_sym1[k, j - 1]
            }
        }
        dw_sym2 = matrix(rnorm(nsim * steps/2, 0, sqrt(dt)), nrow = nsim/2, 
            ncol = steps)
        y_sym2 = matrix(0, nrow = nsim/2, ncol = steps)
        for (j in 1:nrow(y_sym2)) {
            for (l in 2:bk) {
                y_sym2[j, l] <- y_sym2[j, l - 1] + alpha * (theta1 - 
                  y_sym2[j, l - 1]) * dt + sqrt(sig2) * dw_sym2[j, l - 1]
            }
            for (m in (bk + 1):ncol(y_sym2)) {
                y_sym2[j, m] <- y_sym2[j, m - 1] + alpha * (theta2 - 
                  y_sym2[j, m - 1]) * dt + sqrt(sig2) * dw_sym2[j, m - 1]
            }
        }
        # matrices for allopatric pairs
        dw_allo1 = matrix(rnorm(nsim * steps/2, 0, sqrt(dt)), nrow = nsim/2, 
            ncol = steps)
        y_allo1 = matrix(0, nrow = nsim/2, ncol = steps)
        for (k in 1:nrow(y_allo1)) {
            for (i in 2:ncol(y_allo1)) {
                y_allo1[k, i] <- y_allo1[k, i - 1] + alpha * (theta3 - y_allo1[k, 
                  i - 1]) * dt + sqrt(sig2) * dw_allo1[k, i - 1]
            }
        }
        dw_allo2 = matrix(rnorm(nsim * steps/2, 0, sqrt(dt)), nrow = nsim/2, 
            ncol = steps)
        y_allo2 = matrix(0, nrow = nsim/2, ncol = steps)
        for (j in 1:nrow(y_allo2)) {
            for (l in 2:ncol(y_allo2)) {
                y_allo2[j, l] <- y_allo2[j, l - 1] + alpha * (theta1 - 
                  y_allo2[j, l - 1]) * dt + sqrt(sig2) * dw_allo2[j, l - 1]
            }
        }
        y = list(y_sym1, y_sym2, y_allo1, y_allo2)
    }
    if (plot == TRUE) {
        if (model == "BM_null" | model == "OU_null") {
            if (is.null(ylim)) {
            	plot(1:steps, y[1, ], ylim = c(min(y), max(y)),
                  type = "l", ann = FALSE, col=col[1], ...)
            } else {
            	plot(1:steps, y[1, ], ylim = ylim, type= "l",
            	  ann = FALSE, col=col[1], ...)
            }
            for (i in 2:nrow(y)) {
                if (i%%2 != 0) 
                  lines(1:steps, y[i, ], col=col[1])
                if (i%%2 == 0) 
                  lines(1:steps, y[i, ], col =col[2])
            }
        }
        else {
        	if (is.null(ylim)) {
              plot(1:steps, y[[1]][1, ], ylim = c(min(y[[1]], y[[2]]), 
              	max(y[[1]], y[[2]])), type = "l", ann = FALSE, col=col[1], ...)
            } else {
                plot(1:steps, y[[1]][1, ], ylim = ylim, type = "l", 
                  ann = FALSE, col=col[1], ...)
            }
            lines(1:steps, y[[2]][1, ], col = col[2])
            for (i in 2:(nsim/2)) {
              lines(1:steps, y[[1]][i, ], col=col[1])
              lines(1:steps, y[[2]][i, ], col = col[2])
            }
        }
        if (labels == TRUE) {
            title(xlab = "Time Steps")
            title(ylab = ("Phenotype"))
        }
    }
    return(y)
}
