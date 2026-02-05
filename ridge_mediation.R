ridge_mediation <- function(data,
                            exposure,
                            outcome,
                            mediator_cols,           
                            hyperpar,
                            tune = TRUE,
                            confounder_cols = NULL,  
                            n_iter = 30,
                            tol = 1e-6,
                            nfold = 3,
                            seed = NULL,
                            every = 10,
                            corr_sign = "pos",
                            verbose = TRUE) {
  
  stopifnot(all(c(exposure, outcome) %in% names(data)))
  
  if (is.null(mediator_cols)) stop("No mediators provided")
  
  # Ensure mediators do not include X, Y, or confounders
  mediator_cols <- setdiff(mediator_cols, c(exposure, outcome, confounder_cols))
  stopifnot(all(mediator_cols %in% names(data)))
  
  if (!is.null(confounder_cols)) {
    stopifnot(all(confounder_cols %in% names(data)))
  }
  
  # Exposure and outcome
  X <- as.numeric(data[, exposure])
  Y <- as.numeric(data[, outcome])
  
  # Mediator matrix 
  M <- as.matrix(data[, mediator_cols, drop = FALSE])
  
  if (!is.numeric(M)) {
    M <- apply(M, 2, function(col) as.numeric(col))
    M <- as.matrix(M)
  }
  
  # Confounders 
  C <- NULL
  if (!is.null(confounder_cols) && length(confounder_cols) > 0) {
    C <- data[, confounder_cols, drop = FALSE]
  }
  
  n <- nrow(data)
  
  if (!is.null(seed)) set.seed(seed)
  
  my_folds <- caret::createFolds(y = Y, k = nfold, list = TRUE, returnTrain = TRUE)
  
  # Initialize z (latent mediator) 
  
  z <- as.numeric(scale(rnorm(n)))
  df <- data.frame(Y = Y, z = z, X = X)
  if (!is.null(C)) df <- cbind(df, C)
  
  # Trace of parameter estimates across iterations
  trace <- data.frame(
    iter = integer(0),
    alpha = numeric(0),
    beta = numeric(0),
    gamma = numeric(0),
    indirect = numeric(0)
  )
  
  if (!corr_sign %in% c("pos", "neg")) {
    stop("corr_sign must be 'pos' or 'neg'")
  }
  
  # Pre-build formulas
  f1 <- if (is.null(C)) {
    as.formula(z ~ X)
  } else {
    as.formula(paste("z ~ X +", paste(confounder_cols, collapse = " + ")))
  }
  
  f2 <- if (is.null(C)) {
    as.formula(Y ~ z + X)
  } else {
    as.formula(paste("Y ~ z + X +", paste(confounder_cols, collapse = " + ")))
  }
  
  # Block maximization loop 
  for (it in seq_len(n_iter)) {
    
    if (verbose && (it %% every == 0 || it == 1)) {
      cat(sprintf("\n[Iter %d/%d]\n", it, n_iter))
    }
    
    # Enforce desired sign between z and Y
    # (or sign of z coefficient when confounders are present)
    if (is.null(C)) {
      flip <- if (corr_sign == "pos") cor(df$z, df$Y) < 0 else cor(df$z, df$Y) > 0
      if (isTRUE(flip)) df$z <- -df$z
    } else {
      fit2_tmp <- lm(f2, data = df)
      z_coef <- unname(coef(fit2_tmp)["z"])
      flip <- if (corr_sign == "pos") z_coef < 0 else z_coef > 0
      if (isTRUE(flip)) df$z <- -df$z
    }
    
    df$z <- as.numeric(scale(df$z))
    
    # Model 1: z ~ X (+ C)
    fit1 <- lm(f1, data = df)
    coefs1 <- coef(fit1)
    
    # Full linear predictor from model 1
    Xmat1 <- model.matrix(fit1)
    h <- as.numeric(Xmat1 %*% coefs1)
    
    alpha <- unname(coefs1["X"])
    
    # Model 2: Y ~ z + X (+ C)
    fit2 <- lm(f2, data = df)
    coefs2 <- coef(fit2)
    
    beta  <- unname(coefs2["z"])
    gamma <- unname(coefs2["X"])
    
    # Residuals after removing the effect of z
    Xmat2 <- model.matrix(fit2)
    coefs2_noz <- coefs2
    coefs2_noz["z"] <- 0
    linpred_noz <- as.numeric(Xmat2 %*% coefs2_noz)
    e <- Y - linpred_noz
    
    # Target for ridge regression
    d <- (beta * e + h) / (beta^2 + 1)
    
    train_control <- caret::trainControl(
      method = ifelse(isTRUE(tune), "cv", "none"),
      index = my_folds,
      verboseIter = FALSE
    )
    
    ridge_model <- caret::train(
      y = d,
      x = M,
      method = "glmnet",
      trControl = train_control,
      tuneGrid = hyperpar,
      preProcess = c("center", "scale")
    )
    
    best_lambda <- ridge_model$bestTune$lambda
    
    # Update z = Φ(M)
    df$z <- as.numeric(predict(ridge_model, newdata = M))
    
    # Store parameter trace
    trace <- rbind(trace, data.frame(
      iter = it,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      indirect = alpha * beta
    ))
    
    # Convergence check
    if (it > 2) {
      curr <- as.numeric(trace[it, c("alpha", "beta", "gamma", "indirect")])
      prev <- as.numeric(trace[it - 1, c("alpha", "beta", "gamma", "indirect")])
      
      rel_change <- abs(prev - curr) / pmax(abs(prev), .Machine$double.eps)
      
      if (all(rel_change < tol)) {
        cat(sprintf("\n[Iter %d/%d]", it, n_iter))
        cat(sprintf(
          "\nalpha=%.4f  beta=%.4f  gamma=%.5f  (alpha*beta)=%.4f  | bestCV: lambda=%.4f",
          alpha, beta, gamma, alpha * beta, best_lambda
        ))
        break
      }
    }
    
    if (verbose && (it %% every == 0 || it == 1)) {
      cat(sprintf(
        " alpha=%.4f  beta=%.4f  gamma=%.5f  (alpha*beta)=%.4f  | bestCV: lambda=%.4f",
        alpha, beta, gamma, alpha * beta, best_lambda
      ))
    }
  }
  
  list(
    final_model     = ridge_model,   
    final_z         = df$z,           
    path_trace      = trace,          
    mediator_cols   = mediator_cols,
    confounder_cols = confounder_cols,
    fit1            = fit1,
    fit2            = fit2
  )
}
