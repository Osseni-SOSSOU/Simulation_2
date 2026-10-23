## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----warning=FALSE, message=FALSE------------------------------------------
# Chargement des bibliothèques nécessaires
library(fda)          # Analyse fonctionnelle (B-splines et Fourier)
library(compositions) # Analyse de données compositionnelles
library(Compositional)
library(tidyverse)    # Manipulation et visualisation
library(gridExtra)    # Arrangement de graphiques
library(pracma)       # Fonctions mathématiques supplémentaires
library(writexl)    # Pour exporter les donnés en format excel
library(readxl)       # Pour importer les donnée en format excel
library(ggplot2)      # For visualization
library(viridis)      # For color palettes
library(RANN)         # For find k nearest neighbors
library(splines)
library(splines2)


## --------------------------------------------------------------------------

create_bspline_basis <- function(t_grid, n_basis, degree){
  #' @param t_grid observation points of the curve
  #' @param n_bases number of B-splines
  #' @param degree degree of the splines
  
  base <- splines2::bSpline(t_grid, df = n_basis, degree = degree)
  
  return(list(
      base = base,
      n_basis = n_basis
    ))
}


## --------------------------------------------------------------------------
Add_log_ratio <- function(f){
  #' @param f
  #' @return add log-ratio of f
  
  f_alr <- matrix(0, nrow(f)-1, ncol(f))
  for (j in 1:nrow(f_alr)) {
    f_alr[j, ] <- log((f[j+1, ]/f[1,]))
  }
  
  return(f_alr)
}


## --------------------------------------------------------------------------
inverse_alr  <- function(z) {
  #' @param z
  #' @return inverse of z
  
  z_inv <- matrix(0, nrow(z)+1, ncol(z))
  
  z_inv[1, ] <- 1/(1 + colSums(exp(z)))
  
  for (j in 2:(nrow(z)+1)) {
    z_inv[j, ] <- (exp(z[(j-1), ]))/(1 + colSums(exp(z)))
  }
  
  return(z_inv)
}


## --------------------------------------------------------------------------
Kulback_Leibler_divergence <- function(P, Q, t_grid){
  #' @param P
  #' @param Q
  #' @param t_grid
  #' @return Kulback Leibler divergence
  
  KL <- 0
  for (j in 1:nrow(P)) {
    P[P == 0] <- 1e-10
    Q[Q == 0] <- 1e-10
    term <- P[j, ]*log(P[j, ]/Q[j, ])
    # Approximate the integral using the trapezoidal rule.
    integ <- pracma::trapz(t_grid, term)
    
    KL <- KL + integ
  }
  return(KL)
}


## --------------------------------------------------------------------------
Jensen_divergence <- function(P, Q, t_grid, geom = TRUE){
  #' @param P
  #' @param Q
  #' @param t_grid
  #' @return Jensen divergence 
  
  
  if (geom == TRUE){
    delta <- sqrt(P * Q)
  }
  else {
    delta <- (P + Q) / 2
  }
  
  Js <- Kulback_Leibler_divergence(P, delta, t_grid) + Kulback_Leibler_divergence(Q, delta, t_grid)
  
  return(0.5*Js)
}


## --------------------------------------------------------------------------
# Fréchet mean estimator in the functional simplex

frechet_mean <- function(Y_list, alpha) {
  #' @param Y_list List of compositional functions
  #' @param alpha Parameter 
  #' @return The Fréchet mean

  N <- length(Y_list)
  D <- nrow(Y_list[[1]])
  n_points <- ncol(Y_list[[1]])
  
  if (alpha == 0) {
    # Limit case: geometric mean
    log_sum <- matrix(0, nrow = D, ncol = n_points)
    for (i in 1:N) {
      Y_list[[i]][Y_list[[i]] == 0] <- 1e-10
      log_sum <- log_sum + log(Y_list[[i]])
    }
    geometric_mean <- exp(log_sum / N)
    # Closure
    frechet_mean <- geometric_mean / colSums(geometric_mean)
  } 
  else {
    sum_alpha <- matrix(0, nrow = D, ncol = n_points)
    for (i in 1:N) {
      sum_alpha <- sum_alpha + (Y_list[[i]])^alpha
    }
    frechet_mean <- (sum_alpha)^(1/alpha)
    # Closure
    frechet_mean <- frechet_mean / colSums(frechet_mean)
    frechet_mean[is.na(frechet_mean) | is.infinite(frechet_mean)] <- 0
  }
  
  return(frechet_mean)
}


## --------------------------------------------------------------------------
# Estimator of alpha-k-NN

alpha_knn_estimator <- function(X, Y_list, x_new, k, alpha) {
  #' @param X Vector of predictors
  #' @param Y_list List of compositional response functions
  #' @param x_new New predictor value
  #' @param k Number of nearest neighbors
  #' @param alpha Parameter 
  #' @return Estimator 

  
  # 1. Find the k nearest neighbors
  
  if(is.matrix(X)){
    n <- nrow(X)
    distances <- numeric(nrow(X))
    for (i in 1:nrow(X)) {
      distances[i] <- sqrt(sum((X[i, ]-x_new)^2))
    }
  }
  else{
    n <- length(X)
    distances <- abs(X - x_new)
  }
  k_nearest <- order(distances)[1:min(k, n)]
  
  # 2. Extract the corresponding responses
  Y_knn <- Y_list[k_nearest]
  
  # 3. Compute the Fréchet mean
  estimate <- frechet_mean(Y_knn, alpha)
  
  return(estimate)
}


## --------------------------------------------------------------------------
# Estimator of alpha-kernel

# Gaussian kernel
kernel_gaussian <- function(u,h) {
  #' @param u
  #' @param h
  #' @return value
  
  exp(-0.5 * (u/h)^2) / sqrt(2 * pi)
}

# Estimator

alpha_kernel_estimator <- function(X, Y_list, x_new, kernel_funct, h, alpha){
  #' @param X Vector of predictors
  #' @param Y_list List of compositional response functions
  #' @param x_new New predictor value
  #' @param kernel_funct
  #' @param h parameter
  #' @param alpha Parameter 
  #' @return Estimator de m(x_new)
  
  n <- length(Y_list)
  distance <- numeric(n)
  
  # Calcul des distances
  for (i in 1:n) {
    if(is.matrix(X)){
    for (i in 1:nrow(X)) {
      distance[i] <- sqrt(sum((X[i, ]-x_new)^2))
    }
  }
  else{
    distance[i] <- abs(X[i] - x_new)
  }
  }
  
  # Calcul des poids
  weig <- kernel_funct(distance, h)
  
  # Estimation
  Y_est <- matrix(0, nrow = nrow(Y_list[[1]]), ncol = ncol(Y_list[[1]]))
  for (i in 1:n) {
    Y_est <- Y_est + weig[i]*(Y_list[[i]]^(alpha))
  }
  
  # Normalisation
  Y_est <-  Y_est / colSums(Y_est)
  Y_est[is.na(Y_est) | is.infinite(Y_est)] <- 0
  
  return(Y_estimate = Y_est)
}


## --------------------------------------------------------------------------
# KLD regression using multinomial logit
kld_regression_predict <- function(X_train, Y_train, X_test) {
  #'@param X_train vecteur des prédicteurs pour l'entraînement
  #'@param Y_train liste d'observations, chaque élément est une matrice D x n_points
  #'@param X_test vecteur des prédicteurs pour le test
  
  N <- length(Y_train)
  n_points <- ncol(Y_train[[1]])
  D <- nrow(Y_train[[1]])
  
  Y_array <- array(0, dim = c(D, n_points, N))
  for (i in 1:N) {
    
    Y_array[, , i] <- Y_train[[i]]
  }
  
  Y_pred <- matrix(NA, nrow = D, ncol = n_points)
  X <- X_train
  
  for (t in 1:n_points) {
    Y <- t(Y_array[ , t, ])
    model <- kl.compreg(Y, X, B = 1, ncores = 1, xnew = X_test)
    
    Y_pred[ ,t] <- model$est
  }
  
  Y_pred <- as.matrix(Y_pred)
  return(Y_pred)
  }


## --------------------------------------------------------------------------
# Cross-validation to choose alpha and k or h

cross_validation <- function(X, Y_list, alphas = c(0, 0.5, 1), order_paras, n_folds = 10, mesure, t, mode = 0) {
  #' @param X Vector of predictors
  #' @param Y_list List of compositional functions
  #' @param alphas Vector of alpha values to test
  #' @param order_paras
  #' @param n_folds Number of folds
  #' @param mesure
  #' @param t
  #' @param mode
  #' @return The best pair (alpha, k) and the scores.
  
  n <- if(is.matrix(X)) nrow(X) else length(X)
  
  # Create the folds for cross-validation.
  fold_indices <- cut(seq(1, n), breaks = n_folds, labels = FALSE)
  folds <- sample(fold_indices)  # Random shuffle.
  
  # Matrix to store the scores.
  scores <- array(NA, dim = c(length(alphas), length(order_paras), n_folds))
  
  for (fold in 1:n_folds) {
    # Training and test indices.
    test_indices <- which(folds == fold)
    train_indices <- which(folds != fold)
    
    if(is.matrix(X)){
      X_train <- X[train_indices, ]
      X_test <- X[test_indices, ]
      n_test <- nrow(X_test)
    }
    else{
      X_train <- X[train_indices]
      X_test <- X[test_indices]
      n_test <- length(X_test)
    }
    
    Y_train <- Y_list[train_indices]
    Y_test <- Y_list[test_indices]
    
    for (i in seq_along(alphas)) {
      for (j in seq_along(order_paras)) {
          alpha <- alphas[i]
          order_para <- order_paras[j]
        
        
        # Predictions on the test set.
        predictions <- list()
        
          
        if(mode == 0){
          for (idx in 1:n_test) {
          pred <- if(is.matrix(X)) alpha_knn_estimator(X_train, Y_train, X_test[idx,], order_para, alpha) else alpha_knn_estimator(X_train, Y_train, X_test[idx], order_para, alpha)
          
          predictions[[idx]] <- pred
          }
        }
        else{
          for (idx in 1:n_test) {
          pred <- if(is.matrix(X)) alpha_kernel_estimator(X_train, Y_train, X_test[idx,], kernel_gaussian, order_para, alpha) else alpha_kernel_estimator(X_train, Y_train, X_test[idx], kernel_gaussian, order_para, alpha)
          
          predictions[[idx]] <- pred
          }
        }
        
        # Calculation of the error.
        error <- 0
        for (idx in 1:n_test)
          {
          js <- mesure(predictions[[idx]], Y_test[[idx]], t)
          
          error <- error + js
        }
        
        scores[i, j, fold] <- error / length(Y_test)
      }
    }
  }
  
  # Average of the scores over the folds.
  mean_scores <- apply(scores, c(1, 2), mean, na.rm = TRUE)
  
  # Find the best pair (alpha, k).
  best_idx <- which(mean_scores == min(mean_scores), arr.ind = TRUE)
  best_alpha <- alphas[best_idx[1]]
  best_order_para <- order_paras[best_idx[2]]
  best_para <- c(best_alpha, best_order_para)
  
  return(list(
    best_parameter = best_para,
    mean_scores = mean_scores,
    scores = scores,
    alphas = alphas,
    order_paras = order_paras
  ))
}


## --------------------------------------------------------------------------
# Empirical study of divergence

study_divergence <- function(sample_sizes, n_replications, x, alphas, order_paras, D=3, p,  mode = 0, zero_prop = 0, v = 1, seg = TRUE) {
  #' @param sample_sizes
  #' @param n_replications
  #' @param x
  #' @param order_paras
  #' @param D
  #' @param p
  #' @param mode
  #' @return divergence
  
  divergence_size <- matrix(NA, nrow = length(sample_sizes), ncol = n_replications)
  
  #zero_prop <- zero_prop
  #seg <- seg
  
  for (i in seq_along(sample_sizes)) {
    N_current <- sample_sizes[i]
    cat(sprintf(" Sample size: %d\n", N_current))
    
    for (rep in 1:n_replications) {
      # Generate data.
      sim_data <- generate_simulation_sample(N_current, D, p, t_grid, Base, zero_prop, v, seg)
      
      # Cross-validation to save time.
      cv_result <- cross_validation(
              sim_data$X_n,
              sim_data$Y_list,
              alphas = alphas,
              order_paras = order_paras,
              n_folds = 10,
              Kulback_Leibler_divergence,
              t_grid,
              mode
              )
      
      # Estimation at point x.
      true_val <- generat_response_i(x, D, p, t_grid, Base, v, seg)
      
      if(mode == 0){
          estimate <- alpha_knn_estimator(sim_data$X_n, sim_data$Y_list, x, cv_result$best_parameter[2], cv_result$best_parameter[1])
        }
        else{
          estimate <- alpha_kernel_estimator(sim_data$X_n, sim_data$Y_list, x, kernel_gaussian, cv_result$best_parameter[2], cv_result$best_parameter[1])
        }
      
      # Calculation of the error.
      
      mesure_error <- Kulback_Leibler_divergence(estimate, true_val$Y_comp, t_grid)
      
      divergence_size[i, rep] <- mesure_error
    }
  }
  
  # Calculation of the mean errors.
  mean_divergence <- rowMeans(divergence_size)
  
  df_divergence <- data.frame(
    n = sample_sizes,
    divergence = mean_divergence
  )
  
  return(df_divergence)
}


## --------------------------------------------------------------------------
divergence_comp <- function(sample_sizes, n_replications, x, alphas, order_paras, Ds, p, mode = 0, zero_prop = 0, v = 1, seg = TRUE){
  #' @param sample_sizes
  #' @param n_replications
  #' @param x
  #' @param alphas
  #' @param order_paras
  #' @param Ds
  #' @param mode
  #' @return liste of divergence
  
    data_divergence <- data.frame(size = sample_sizes)
  number_comp <- c("size")
  
  #zero_prop <- zero_prop
  #seg <- seg
  
  cat("\n Divergence study (this may take a few minutes)... \n")
  for (l in 1:length(Ds)) {
    number_comp[l+1] <- sprintf("D = %d", Ds[l])
    D <- Ds[l]
    
    cat(sprintf("\n Pour: %d composantes \n \n", D))
    P <- study_divergence(sample_sizes, n_replications, x, alphas, order_paras, D, p, mode = mode, zero_prop, v, seg)
    
    data_divergence[ ,l+1] <- P[,2]
    names(data_divergence) <- number_comp
  }
  return(data_divergence)
}


## --------------------------------------------------------------------------
divergence_kld <- function(sample_sizes, n_replications, x, D, p, zero_prop = 0, v = 1, seg = TRUE){
  
  
  divergence_size <- matrix(NA, nrow = length(sample_sizes), ncol = n_replications)
  
  zero_prop <- zero_prop
  seg <- seg
  
  for (i in seq_along(sample_sizes)) {
    N_current <- sample_sizes[i]
    cat(sprintf(" Sample size: %d\n", N_current))
    
    for (rep in 1:n_replications) {
      sim_data <- generate_simulation_sample(N_current, D, p, t_grid, Base, zero_prop, v, seg)
      
      true_val <- generat_response_i(x, D, p, t_grid, Base, v, seg)
      
       # Estimation at point x.
      estimate <- kld_regression_predict(sim_data$X_n, sim_data$Y_list, x)
      
      # Calculation of the error.
      
      mesure_error <- Kulback_Leibler_divergence(estimate, true_val$Y_comp, t_grid)
      
      divergence_size[i, rep] <- mesure_error
      
    }
  }
  
  # Calculation of the mean errors.
  mean_divergence <- rowMeans(divergence_size)
  
  df_divergence <- data.frame(
    n = sample_sizes,
    divergence = mean_divergence
  )
  
  return(df_divergence)
}


## --------------------------------------------------------------------------
divergence_kld_comp <- function(sample_sizes, n_replications, x, Ds, p, zero_prop = 0, v = 1, seg = TRUE){
  
  data_divergence <- data.frame(size = sample_sizes)
  number_comp <- c("size")
  
  zero_prop <- zero_prop
  seg <- seg
  
  cat("\n Divergence study (this may take a few minutes)... \n")
  for (l in 1:length(Ds)) {
    number_comp[l+1] <- sprintf("D = %d", Ds[l])
    D <- Ds[l]
    
    cat(sprintf("\n Pour: %d composantes \n \n", D))
    P <- divergence_kld(sample_sizes, n_replications, x, D, p, zero_prop, v, seg)
    
    data_divergence[ ,l+1] <- P[,2]
    names(data_divergence) <- number_comp
  }
  return(data_divergence)
}


## --------------------------------------------------------------------------
# Simulation de la réponse i

generat_X <- function(p, n = 1, seg = TRUE){
  #' @param n
  #' @param p
  
  if(seg){
    X <- numeric(0)
    for (i in 1:n) {
      X[i] <- runif(1, min = -1, max = 1)
    }
  }
  else{
    X <- matrix(NA, nrow = n, ncol = p)
    for (i in 1:n) {
    X[i, ] <- rnorm(p)
    }
  }
  
  return(X)
}

generat_response_i <- function(x, D, p, t_grid, basis_info, v = 1, seg = TRUE){
  #' @param x vecteur des prédicteurs
  #' @param t_grid grille temporelle
  #' @param p
  #' @param v
  #' @param seg
  #' @return liste contenant les données générées
  
  n_basis <- basis_info$n_basis
  n_grid <- length(t_grid)
  
  
  y_i <- matrix(0, nrow = D-1, ncol = n_grid)
  for (j in 1:(D-1)) {
    
    # Génération d'une erreur
      noise_coef <- rnorm(n_basis, 0, 0.05)
      epsilon_smooth <- basis_info$base %*% noise_coef
    
    # Application du modèle
      
      if(seg){
        #Génération des coefficients
        
        beta_j <- matrix(NA, nrow = 2, ncol = n_grid)
        
          beta_coef_1 <- rnorm(n_basis, -1, 0.3)
          beta_coef_2 <- rnorm(n_basis, 1, 0.2)
          
          beta_j[1, ] <- basis_info$base %*% beta_coef_1
          beta_j[2, ] <- basis_info$base %*% beta_coef_2
        
        if(x > 0){
          y_i[j, ] <- (x^2) * beta_j[1, ] + epsilon_smooth
        }
        else{
          y_i[j, ] <- (x^3) * beta_j[2, ] + epsilon_smooth
        }
      }
      else{
        #Génération des coefficients
        beta_0_coef <- rnorm(n_basis, 2, 0.5)
        beta_0 <- basis_info$base %*% beta_0_coef
      
        beta_j <- matrix(NA, nrow = p, ncol = n_grid)
        for (r in 1:p) {
          beta_coef <- rnorm(n_basis, -3, 1)
          beta_j[r, ] <- basis_info$base %*% beta_coef
        }
        
        y_i[j, ] <- beta_0 + as.vector((x^v) %*% beta_j) + epsilon_smooth
      }
      
  }
  
    # Retour sur le simplexe fonctionnelle
    y_comp <- inverse_alr(y_i)
  
    return(list(
      Y_fonc = y_i,
      Y_comp = y_comp
    ))
    
}


## --------------------------------------------------------------------------
generate_simulation_sample <- function(n, D, p, t_grid, basis_info, zero_prop = 0, v = 1, seg = TRUE){
  #' @param n
  #' @param X
  #' @param t_grid
  #' @param basis_info
  #' @return liste contenant les données générées
  
  
  Y_n <- list()
  X_n <- generat_X(p, n, seg)
  
  seg <- seg
  
  # Introduce zeros if specified (at random time points and components)
  if (zero_prop > 0) {
    Y_array <- array(0, dim = c(D, length(t_grid), n))
    
    for (i in 1:n) {
      if(is.matrix(X_n)){
        y_i <- generat_response_i(X_n[i, ], D, p, t_grid, basis_info, v, seg)
        Y_array[ , , i] <- y_i$Y_comp
      }
      else{
        y_i <- generat_response_i(X_n[i], D, p, t_grid, basis_info, v, seg)
        Y_array[ , , i] <- y_i$Y_comp
      }
    }
    
    n_zeros <- round(n * length(t_grid) * D * zero_prop)
    zero_indices <- sample(1:(n * length(t_grid) * D), n_zeros)
    Y_array[zero_indices] <- 0
    
    # Re-normalize each observation at each time point
    for (i in 1:n) {
      for (t in 1:length(t_grid)) {
        if (sum(Y_array[, t, i]) > 0) {
          Y_array[, t, i] <- Y_array[, t, i] / sum(Y_array[, t, i])
          }
        else {
          Y_array[, t, i] <- rep(1/D, D)
          }
      }
      Y_n[[i]] <- Y_array[, , i]
    }
  }
  else{
   for (i in 1:n) {
     if(is.matrix(X_n)){
        y_i <- generat_response_i(X_n[i, ], D, p, t_grid, basis_info, v, seg)
        Y_n[[i]] <- y_i$Y_comp
      }
      else{
        y_i <- generat_response_i(X_n[i], D, p, t_grid, basis_info, v, seg)
        Y_n[[i]] <- y_i$Y_comp
      }
    } 
  }
  
  return(list(
    Y_list = Y_n,
    X_n = X_n
  ))
}


## --------------------------------------------------------------------------
# Visualization of a typical compositional function.
plot_compositional_function <- function(Y, t, title = "Compositional function.") {
  #' @param Y
  #' @param t
  #' @param title
  
  D <- nrow(Y)
  df <- data.frame(
    t = rep(t, D),
    value = as.vector(t(Y)),
    component = factor(rep(1:D, each = length(t)))
  )
  
  a <- min(Y) + 0.02
  b <- max(Y) + 0.02
  
  p <- ggplot(df, aes(x = t, y = value, color = component)) +
    geom_line(linewidth = 1) +
    theme_minimal() +
    labs(
      title = title,
      x = "Time (t)",
      y = "Proportion",
      color = "Composante"
    ) +
    scale_color_viridis(discrete = TRUE) +
    ylim(a, b) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right"
    )
  
  return(p)
}


## --------------------------------------------------------------------------
# Visualization of divergence.
plot_divergence_function <- function(data_div, title = "function") {
  #' @param Y
  #' @param title
  
  a <- min(data_div[,-1]) - 0.03
  b <- max(data_div[,-1]) + 0.03
  
  t <- data_div[,1]
  Y_lim = ylim(a, b)
  
  df <- data.frame(
    t = rep(data_div[,1], ncol(data_div)-1),
    value = as.vector(as.matrix(data_div[,-1])),
    number_component = factor(rep(names(data_div)[-1], each = length(data_div[,1])), levels = names(data_div)[-1])
  )
  
  p <- ggplot(df, aes(x = t, y = value, color = number_component)) +
    geom_line(linewidth = 0.7) +
    theme_minimal() +
    labs(
      title = title,
      x = "size",
      y = "Ratio of divergence",
      color = "Nb composante"
    ) +
    scale_color_viridis(discrete = TRUE) +
    Y_lim +
    theme(
      plot.title = element_text(hjust = 1, size = 14),
      legend.position = "right"
    )
  
  return(p)
}


## --------------------------------------------------------------------------
divergnece_ratio <- function(dive1, dive2){
  
  report <- data.frame(size = dive1[,1])
  nr <- nrow(dive1)
  nc <- ncol(dive1)
  
  for (i in 1:nr) {
    for(j in 2:nc) {
      report[i, j] <- dive1[i,j] / dive2[i,j]
    }
  }
  names(report) <- names(dive1)
  return(report)
}

