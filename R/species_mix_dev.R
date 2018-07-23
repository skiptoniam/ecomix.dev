#' @useDynLib ecomix.dev
# #' @importFrom Rcpp sourceCpp
NULL

## need to put disps in the starting vals even if NULL.

sam_optimise <- function(y, X, offset, spp_wts, site_spp_wts, y_is_na, nS, nG, nObs, disty, start_vals, control) {

  inits <- c(start_vals$alphas, start_vals$betas, start_vals$pis, start_vals$disp)
  np <- as.integer(ncol(X[,-1]))
  n <- as.integer(nrow(X))

  # parameters to optimise
  alpha <- as.numeric(start_vals$alphas);
  beta <- as.numeric(start_vals$betas);
  eta <- as.numeric(additive_logistic(start_vals$pis,TRUE)[seq_len(nG-1)])
  disp <- as.numeric(start_vals$disp)

  #scores
  alpha.score <- as.numeric(rep(NA, length(alpha)))
  beta.score <- as.numeric(rep(NA, length(beta)))
  eta.score <- as.numeric(rep(NA, length(eta)))
  disp.score <- as.numeric(rep(NA, length(disp)))
  getscores <- 1
  scores <- as.numeric(rep(NA,length(c(alpha,beta,eta,disp))))

  #model quantities
  pis_out <- as.numeric(rep(NA, nG))  #container for the fitted RCP model
  mus <- as.numeric(array( NA, dim=c( nObs, nS, nG)))  #container for the fitted spp model
  loglikeS <- as.numeric(rep(NA, S))
  loglikeSG  <- as.numeric(matrix(NA, nrow = nS, ncol = nG))

  #c++ call to optimise the model (needs pretty good starting values)
  tmp <- .Call("species_mix_cpp",
               as.numeric(as.matrix(y)), as.numeric(as.matrix(X[,-1])), as.numeric(offset), as.numeric(spp_wts),
               as.numeric(as.matrix(site_spp_wts)), as.integer(!y_is_na),
               # SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rspp_wts, SEXP Rsite_spp_wts, SEXP Ry_not_na, // data
               as.integer(nS), as.integer(nG), as.integer(np), as.integer(nObs), as.integer(disty),
               # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
               as.double(alpha), as.double(beta), as.double(eta), as.double(disp),
               # SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
               alpha.score, beta.score, eta.score, disp.score, as.integer(control$getscores), scores,
               # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
               pis_out, mus, loglikeS, loglikeSG,
               # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
               as.integer(control$maxit_cpp), as.integer(control$trace_cpp), as.integer(control$nreport_cpp),
               as.numeric(control$abstol_cpp), as.numeric(control$reltol_cpp), as.integer(control$conv_cpp), as.integer(control$printparams_cpp),
               # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
               as.integer( control$optimise_cpp), as.integer(control$loglOnly_cpp), as.integer( control$derivOnly_cpp), as.integer(control$optiDisp),
               # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
               PACKAGE = "ecomix.dev")

  ret <- tmp
  ret$mus <- array(mus, dim=c(n, S, G))
  ret$scores <- list(alpha.scores = alpha.score, beta.scores = beta.score, eta.scores=eta.score, disp.scores=disp.score)
  ret$S <- S; ret$G <- G; ret$np <- np; ret$n <- n;
  ret$start.vals <- inits
  ret$loglikeSG <- loglikeSG  #for residuals
  ret$loglikeS <- loglikeS  #for residuals
  return(ret)
}

"additive_logistic" <- function (x,inv=FALSE) {
  if(inv){
    x <- log(x/x[length(x)])
    return(x)
  } else {

    x.t <- exp(x)
    x.t <- x.t/(1+sum(x.t))
    x.t[length(x.t)+1] <- 1-sum(x.t)
    return(x.t)
  }
}
