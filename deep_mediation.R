deep_mediation <- function(data,
                           exposure,
                           outcome,
                           mediator_cols,              # mediator columns
                           hyperpar,                   # data.frame with layer1/layer2/layer3/decay/lr/maxit
                           allowParallel = FALSE,
                           confounder_cols = NULL,     # e.g. c("age", "sex")
                           n_iter = 30,
                           tol = 1e-6,
                           nfold = 3,
                           every = 10,
                           tune = TRUE,
                           different_start = TRUE,
                           corr_sign = "pos",
                           nthreads = NULL,
                           seed = NULL,
                           maxit = 50,
                           stop_tune = NULL,
                           verbose = TRUE) {
  
  # Custom caret model: RSNNS MLP with tunable learning rate (lr) and decay
  mlpWeightDecayML_custom <- list(
    label = "Multi-Layer Perceptron, multiple layers",
    library = "RSNNS",
    loop = NULL,
    type = c("Regression", "Classification"),
    
    # Parameters: add lr
    parameters = data.frame(
      parameter = c("layer1", "layer2", "layer3", "decay", "lr"),
      class = c("numeric", "numeric", "numeric", "numeric", "numeric"),
      label = c("#Hidden Units layer1", "#Hidden Units layer2",
                "#Hidden Units layer3", "Weight Decay", "Learning Rate")
    ),
    
    grid = function(x, y, len = NULL, search = "grid") {
      if (search == "grid") {
        expand.grid(
          layer1 = ((1:len) * 2) - 1,
          layer2 = 0,
          layer3 = 0,
          decay  = c(0, 10^seq(-1, -4, length = len - 1)),
          lr     = c(0.1, 0.2, 0.3)
        )
      } else {
        data.frame(
          layer1 = sample(2:20, replace = TRUE, size = len),
          layer2 = sample(c(0, 2:20), replace = TRUE, size = len),
          layer3 = sample(c(0, 2:20), replace = TRUE, size = len),
          decay  = 10^runif(len, min = -5, max = 1),
          lr     = runif(len, min = 0.01, max = 0.5)
        )
      }
    },
    
    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
      theDots <- list(...)
      
      # Prevent overriding of these via ...
      theDots <- theDots[!(names(theDots) %in% c("size", "linOut"))]
      
      if ("learnFunc" %in% names(theDots)) {
        theDots$learnFunc <- NULL
        warning("Cannot over-ride 'learnFunc' argument for this model. BackpropWeightDecay is used.")
      }
      
      if ("learnFuncParams" %in% names(theDots)) {
        prms <- theDots$learnFuncParams
        prms[1] <- param$lr
        prms[2] <- param$decay
        warning("Over-riding learning rate and weight decay values in the 'learnFuncParams' argument you passed in. Other values are retained")
      } else {
        prms <- c(param$lr, param$decay, 0, 0)
      }
      
      if (is.factor(y)) {
        y <- RSNNS:::decodeClassLabels(y)
        lin <- FALSE
      } else {
        lin <- TRUE
      }
      
      nodes <- c(param$layer1, param$layer2, param$layer3)
      if (any(nodes == 0)) {
        nodes <- nodes[nodes > 0]
        warning(
          "At least one layer had zero units and were removed. The new structure is ",
          paste0(nodes, collapse = "->"),
          call. = FALSE
        )
      }
      
      args <- list(
        x = x,
        y = y,
        learnFunc = "BackpropWeightDecay",
        learnFuncParams = prms,
        size = nodes,
        linOut = lin
      )
      
      args <- c(args, theDots)
      do.call(RSNNS::mlp, args)
    },
    
    predict = function(modelFit, newdata, submodels = NULL) {
      out <- predict(modelFit, newdata)
      if (modelFit$problemType == "Classification") {
        out <- modelFit$obsLevels[apply(out, 1, which.max)]
      } else {
        out <- out[, 1]
      }
      out
    },
    
    prob = function(modelFit, newdata, submodels = NULL) {
      out <- predict(modelFit, newdata)
      colnames(out) <- modelFit$obsLevels
      out
    },
    
    levels = function(x) x$obsLevels,
    
    tags = c("Neural Network", "L2 Regularization"),
    
    sort = function(x) {
      x[order(x$layer1, x$layer2, x$layer3, -x$decay, -x$lr), ]
    }
  )
  
  # ---- 1) Extract X, Y, M, C ----
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
  my_folds <- caret::createFolds(y = Y, k = nfold, list = TRUE, returnTrain = TRUE)
  
  # ---- 2) Parallel backend  ----
  cl <- NULL
  
  on.exit({
    # Return to sequential backend
    try(foreach::registerDoSEQ(), silent = TRUE)
    # Stop cluster if created
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
  
  # ---- 3) Initialize z and working df ----
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
  
  
  # Set caret seeds
  n_resamples <- length(my_folds)
  seeds_list <- vector("list", n_resamples + 1)
  
  for (i in seq_len(n_resamples)) {
    if (isTRUE(different_start)) {
      seeds_list[[i]] <- stats::rpois(nrow(hyperpar), 1000)
    } else {
      seeds_list[[i]] <- rep(123, nrow(hyperpar))
    }
  }
  seeds_list[[n_resamples + 1]] <- if (isTRUE(different_start)) stats::rpois(1, 1000) else 123
  
  
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
  
  # ---- Train control ----
  train_control <- caret::trainControl(
    method = ifelse(isTRUE(tune), "cv", "none"),
    index = my_folds,
    verboseIter = FALSE,
    seeds = if (isTRUE(tune)) seeds_list else NA,
    allowParallel = use_parallel
  )
  
  best_params0 <- NULL
  same_tune <- 0
  
  # ---- 4) Block maximization loop ----
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
    
    # Enforce desired sign between z and Y (or z coefficient with confounders)
    if (is.null(C)) {
      flip <- if (corr_sign == "pos") cor(df$z, df$Y) < 0 else cor(df$z, df$Y) > 0
      if (isTRUE(flip)) df$z <- -df$z
    } else {
      fit2_tmp <- stats::lm(f2, data = df)
      z_coef <- unname(stats::coef(fit2_tmp)["z"])
      flip <- if (corr_sign == "pos") z_coef < 0 else z_coef > 0
      if (isTRUE(flip)) df$z <- -df$z
    }
    df$z <- as.numeric(scale(df$z))
    
    # Model 1: z ~ X (+ C)
    fit1 <- stats::lm(f1, data = df)
    coefs1 <- stats::coef(fit1)
    
    Xmat1 <- stats::model.matrix(fit1)
    h <- as.numeric(Xmat1 %*% coefs1)
    alpha <- unname(coefs1["X"])
    
    # Model 2: Y ~ z + X (+ C)
    fit2 <- stats::lm(f2, data = df)
    coefs2 <- stats::coef(fit2)
    
    beta  <- unname(coefs2["z"])
    gamma <- unname(coefs2["X"])
    
    Xmat2 <- stats::model.matrix(fit2)
    coefs2_noz <- coefs2
    coefs2_noz["z"] <- 0
    linpred_noz <- as.numeric(Xmat2 %*% coefs2_noz)
    e <- Y - linpred_noz
    
    # Target for the deep model
    d <- (beta * e + h) / (beta^2 + 1)
    
    deep_model <- caret::train(
      x = M,
      y = d,
      method = mlpWeightDecayML_custom,
      trControl = train_control,
      tuneGrid = hyperpar,
      linOut = TRUE,
      maxit = maxit,
      preProcess = c("center", "scale")
    )
    
    best_params <- deep_model$bestTune
    
    if (it > 1 && all(best_params == best_params0) && same_tune < 10) {
      same_tune <- same_tune + 1
    } else if (same_tune < 10) {
      same_tune <- 0
    }
    best_params0 <- best_params
    
    # Update z = Φ(M)
    df$z <- as.numeric(predict(deep_model, newdata = M))
    
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
            "\n[Iter %d/%d]\n alpha=%.4f  beta=%.4f  gamma=%.4f  (alpha*beta)=%.4f  | best CV: decay=%s",
            it, n_iter, alpha, beta, gamma, alpha * beta,
            as.character(best_params$decay)))
        }
        break
      }
    }
    
    if (verbose && (it %% every == 0 || it == 1)) {
      cat(sprintf(
        " alpha=%.4f  beta=%.4f  gamma=%.4f  (alpha*beta)=%.4f  | best CV: decay=%s",
        alpha, beta, gamma, alpha * beta,
        as.character(best_params$decay)))
    }
  }
  
  # ---- Output ----
  list(
    final_model     = deep_model,
    final_z         = df$z,
    path_trace      = trace,
    mediator_cols   = mediator_cols,
    confounder_cols = confounder_cols,
    fit1            = fit1,
    fit2            = fit2
  )
}