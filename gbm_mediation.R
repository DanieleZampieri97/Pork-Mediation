gbm_mediation <- function(data,
                          exposure,
                          outcome,
                          mediator_cols,
                          hyperpar,
                          allowParallel = FALSE,
                          nthreads = NULL,
                          confounder_cols = NULL,
                          n_iter = 30,
                          tol = 1e-6,
                          nfold = 3,
                          seed = NULL,
                          every = 10,
                          tune = TRUE,
                          corr_sign = "pos",
                          verbose = TRUE) {
  
  # Extract X, Y, M, C 
  stopifnot(all(c(exposure, outcome) %in% names(data)))
  if (is.null(mediator_cols)) stop("No mediators provided")
  
  mediator_cols <- setdiff(mediator_cols, c(exposure, outcome, confounder_cols))
  stopifnot(all(mediator_cols %in% names(data)))
  
  if (!is.null(confounder_cols)) {
    stopifnot(all(confounder_cols %in% names(data)))
  }
  
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
  
  n <- nrow(M)
  
  if (!is.null(seed)) set.seed(seed)
  
  # Parallel backend 
  cl <- NULL
  
  on.exit({
    try(foreach::registerDoSEQ(), silent = TRUE)
    if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  use_parallel <- isTRUE(allowParallel) && isTRUE(tune)
  
  if (isTRUE(use_parallel)) {
    if (is.null(nthreads)) {
      max_cores <- parallel::detectCores()
      nthreads <- min(nrow(hyperpar) * nfold, max(1, max_cores - 1))
    }
    cl <- parallel::makeCluster(nthreads)
    doParallel::registerDoParallel(cl)
  }
  
  my_folds <- caret::createFolds(y = Y, k = nfold, list = TRUE, returnTrain = TRUE)
  
  # Initialize z (latent mediator) and df
  z <- as.numeric(scale(rnorm(n)))
  df <- data.frame(Y = Y, z = z, X = X)
  if (!is.null(C)) df <- cbind(df, C)
  
  trace <- data.frame(
    iter = integer(0),
    alpha = numeric(0),
    beta = numeric(0),
    gamma = numeric(0),
    indirect = numeric(0)
  )
  
  # Pre-build formulas
  f1 <- if (is.null(C)) {
    stats::as.formula(z ~ X)
  } else {
    stats::as.formula(paste("z ~ X +", paste(confounder_cols, collapse = " + ")))
  }
  
  f2 <- if (is.null(C)) {
    stats::as.formula(Y ~ z + X)
  } else {
    stats::as.formula(paste("Y ~ z + X +", paste(confounder_cols, collapse = " + ")))
  }
  
  # Train control
  train_control <- caret::trainControl(
    method = ifelse(isTRUE(tune), "cv", "none"),
    index = my_folds,
    verboseIter = FALSE,
    allowParallel = use_parallel
  )
  
  
  best_params0 <- NULL
  same_tune <- 0
  
  # Block maximization loop
  for (it in seq_len(n_iter)) {
    
    if (verbose && (it %% every == 0 || it == 1)) {
      cat(sprintf("\n[Iter %d/%d]\n", it, n_iter))
    }
    
    # If tuning results repeat for 10 iterations, stop tuning
    if (same_tune == 10) {
      hyperpar <- best_params
      tune <- FALSE
      
      train_control <- caret::trainControl(
        method = "none",
        verboseIter = FALSE,
        allowParallel = FALSE
      )
      
      if (isTRUE(allowParallel)) doParallel::registerDoSEQ()
      same_tune <- Inf
    }
    
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
    
    # Residuals after removing the contribution of z
    Xmat2 <- model.matrix(fit2)
    coefs2_noz <- coefs2
    coefs2_noz["z"] <- 0
    linpred_noz <- as.numeric(Xmat2 %*% coefs2_noz)
    e <- Y - linpred_noz
    
    # Target for GBM
    d <- (beta * e + h) / (beta^2 + 1)
    
    gbm_model <- caret::train(
      x = M,
      y = d,
      method = "gbm",
      trControl = train_control,
      tuneGrid = hyperpar,
      metric = "RMSE",
      verbose = FALSE
    )
    
    best_params <- gbm_model$bestTune
    
    if (it > 1 && all(best_params == best_params0) && same_tune < 10) {
      same_tune <- same_tune + 1
    } else if (same_tune < 10) {
      same_tune <- 0
    }
    best_params0 <- best_params
    
    # Update z = Φ(M)
    df$z <- as.numeric(predict(gbm_model, newdata = M))
    
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
        if (verbose) {
          cat(sprintf(
            "  alpha=%.5f  beta=%.5f  gamma=%.6f  (alpha*beta)=%.6f  | best CV: ntrees=%s, inter_depth=%s, shrink=%.3f",
            alpha, beta, gamma, alpha * beta,
            best_params$n.trees,
            best_params$interaction.depth,
            best_params$shrinkage
          ))
        }
        break
      }
    }
    
    if (verbose && (it %% every == 0 || it == 1)) {
      cat(sprintf(
        "  alpha=%.5f  beta=%.5f  gamma=%.6f  (alpha*beta)=%.6f  | best CV: ntrees=%s, inter_depth=%s, shrink=%.3f",
        alpha, beta, gamma, alpha * beta,
        best_params$n.trees,
        best_params$interaction.depth,
        best_params$shrinkage
      ))
    }
  }
  
  # Output 
  list(
    final_model     = gbm_model,
    final_z         = df$z,
    path_trace      = trace,
    mediator_cols   = mediator_cols,
    confounder_cols = confounder_cols,
    fit1            = fit1,
    fit2            = fit2
  )
}