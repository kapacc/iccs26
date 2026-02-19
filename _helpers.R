non_scar_labelling.mvc <- function(df, target_c_calc, n_vars = 1) {
  
  # Remove columns Y and S
  df_filtered <- df[, !names(df) %in% c("Y", "S")]
  
  # Check if n_vars does not exceed the number of columns in df_filtered
  if (n_vars >= ncol(df_filtered)) {
    print(paste("n_vars exceeds the number of available columns. Using", ncol(df_filtered) - 1, "variables instead."))
    n_vars <- ncol(df_filtered) - 1
  }
  
  # Calculate variance for each column
  variances <- sapply(df_filtered, var)
  
  # Sort variances in descending order
  sorted_variances <- sort(variances, decreasing = TRUE)
  
  # Find the names of columns with the highest variances (number depends on n_vars)
  top_variable_columns <- names(sorted_variances)[1:n_vars]
  
  # Create a temporary dataframe that contains the sum of values for selected columns
  df_temp <- df %>%
    mutate(rn = row_number(rowSums(across(all_of(top_variable_columns))))) %>%
    mutate(
      S02 = ifelse(rn / nrow(df) < 0.2, 1, 0) * Y,
      S04 = ifelse(rn / nrow(df) < 0.4, 1, 0) * Y,
      S06 = ifelse(rn / nrow(df) < 0.6, 1, 0) * Y,
      S08 = ifelse(rn / nrow(df) < 0.8, 1, 0) * Y
    )
    
    # Calculate c_calc for three levels
    c_calc02 <- sum(df_temp$S02)/sum(df_temp$Y)
    c_calc04 <- sum(df_temp$S04)/sum(df_temp$Y)
    c_calc06 <- sum(df_temp$S06)/sum(df_temp$Y)
    c_calc08 <- sum(df_temp$S08)/sum(df_temp$Y)
    
    c_calc <- c(c_calc02, c_calc04, c_calc06, c_calc08)
    
    # Fit a linear model for c_calc
    coef <- coef(glm(c_calc ~ c(0.2, 0.4, 0.6, 0.8))) %>% as.numeric()
    
    # Calculate the proportion for target_c_calc
    rn_frac <- (target_c_calc - coef[1])/coef[2]
    
    # Create the final dataframe with the column S
    df1 <- df %>%
      mutate(rn = row_number(rowSums(across(all_of(top_variable_columns))))) %>%
      mutate(S = ifelse(rn / nrow(df) < rn_frac, 1, 0) * Y) %>%
      select(-rn)
      
      # Return the dataframe with the new column S
      return(df1)
}

scar_labelling <- function(df, target_c_calc){
  
  current_frame <- df
  current_frame$Y <- as.numeric(current_frame$Y)
  current_frame$S <- current_frame$Y * rbinom(nrow(current_frame),1,target_c_calc)
  current_frame <- as.data.frame(sapply(current_frame, as.numeric))
  
  return(current_frame)
}

create_correlation_map <- function(data_list, y_var, threshold, output_file) {
  pdf(output_file, width = 10, height = 8)  # Open a PDF file
  
  for (i in seq_along(data_list)) {
    data <- data_list[[i]]
    data_name <- names(data_list)[i]
    
    # Remove the Y column
    data <- data[, !names(data) %in% y_var]
    
    # Select only numeric columns
    numeric_data <- data[, sapply(data, is.numeric), drop = FALSE]
    
    # Check if numeric_data is not empty and if there are at least two numeric columns
    if (is.null(numeric_data)) {
      next
    }
    
    if (ncol(numeric_data) < 2) {
      next
    }
    
    if (ncol(numeric_data) > 1000) {
      next
    }
    
    # Calculate the correlation matrix
    correlation_matrix <- cor(numeric_data, use = "complete.obs")
    
    # Transform the correlation matrix into a long format
    melted_correlation_matrix <- melt(correlation_matrix)
    
    # Add a column to indicate high correlation
    melted_correlation_matrix$high_correlation <- abs(melted_correlation_matrix$value) > threshold
    
    # Create the correlation map
    p <- ggplot(data = melted_correlation_matrix, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = ifelse(high_correlation, "X", "")), color = "black", size = 5) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-1, 1), space = "Lab", 
                           name = "Pearson Correlation Coefficient") +
      theme_minimal() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      coord_fixed() +
      ggtitle(data_name)
    
    print(p)  # Save the plot on a new page in the PDF file
  }
  
  dev.off()  # Close the PDF file
}

find_glm_coef <- function(train_df, pecked_part = 0.25, n_of_pecking = 5, lasso = F, strict = F, keep_original_1 = TRUE) {
  
  for (k in 1:n_of_pecking) {
    set.seed(k)
    q <- pecked_part
    
    S1_ind <- which(train_df$S == 1)
    S0_ind <- which(train_df$S == 0)
    
    S1_ind_sample <- sample(S1_ind, length(S1_ind) * q)
    
    # The set S=1 without q% observations is denoted as C1(-q)
    C1_minus_q <- train_df[setdiff(S1_ind,S1_ind_sample),] 
    
    # The set S=0 with added q% observations from class S=1 is denoted as C0(+q)
    C0_plus_q <- train_df[c(S0_ind,S1_ind_sample),]
    
    # Apply the 2-means algorithm to the set C0(+q).
    C0_plus_q <-
      C0_plus_q %>% 
      mutate_if( ~ length(unique(.)) == 2, as.factor) %>%  # Columns in the dataframe C0_plus_q that have exactly two unique values will be converted to factor type.
      mutate(S = as.numeric(S) - 1, Y = as.numeric(Y) - 1)
    
    # Example of k-means
    mixedClusters <- kmeans(select(C0_plus_q, -c("Y", "S")), centers = 2, nstart = 5)
    table(mixedClusters$cluster)
    C0_plus_q <- as.data.frame(C0_plus_q)
    C0_plus_q$cluster <- as.numeric(mixedClusters$cluster)
    
    # 2. The cluster most similar to C1(-q) obtained from 2-means is denoted as C1 (where most S = 1) and the least similar as C0
    like1clust <-
      C0_plus_q %>%
      group_by(cluster) %>%
      summarise(frac = mean(S)) %>%
      arrange(desc(frac)) %>%
      slice(1) %>% 
      pull(cluster)
    
    C1 <-
      C0_plus_q %>% filter(cluster == like1clust) %>% select(-cluster)
    C0 <-
      C0_plus_q %>% filter(cluster != like1clust) %>% select(-cluster)
    
    # 3. Fit a logistic regression model for Y=1 (prediction of Y) from class C1(-q) and class C1 (from point 2)
    # and for class Y=0 (prediction of Y) for observations from C0 from point 2.
    if(keep_original_1 == FALSE) {
      to_Y_modelling <-
        rbind.data.frame(cbind.data.frame(rbind.data.frame(C1_minus_q, C1), YY = 1),
                         cbind.data.frame(C0, YY = 0)) %>% select(-c(Y, S)) %>% mutate_all(as.numeric)
    } else {
      to_Y_modelling <-
        rbind.data.frame(
          # Our improvement label contains as follows: our unpecked S == 1, our whole 'like1' cluster, our S == 1 returning from 'like0' cluster
          # To sum up: we don't have the intention of losing any of our S==1 labels 
          cbind.data.frame(rbind.data.frame(C1_minus_q, C1, C0 %>% filter (S == 1)), YY = 1),
          cbind.data.frame(C0 %>% filter (S == 0), YY = 0)) %>% select(-c(Y, S)) %>% mutate_all(as.numeric)
    }
    
    if(lasso == T) {
      x_to_lasso <- to_Y_modelling %>% select(-YY) %>% as.matrix
      obj3 <- glmnet::cv.glmnet(x = x_to_lasso, y = to_Y_modelling$YY,standardize=TRUE, intercept=TRUE,family="binomial")
      lambda <- obj3$lambda.min
      delta = lambda * 0.5
      betasx1<-coefficients(obj3,s=lambda)
      lasso_betas <- as.vector(betasx1)
      names(lasso_betas) <- row.names(betasx1) 
      current_coef <- ifelse(abs(lasso_betas) > delta,lasso_betas,0)
    } else {
      glm_clust <- my_glm(YY ~ ., data = to_Y_modelling, family = binomial)
      current_coef <- coef(glm_clust)
    }
    
    if (k == 1) {
      coef_df <- data.frame(row.names = names(current_coef))
      coef_df[[k]] <- current_coef
    } else {
      coef_df[[k]] <- current_coef
    }
  }
  
  if (strict == F) {
    glm_coef <- apply(coef_df, 1, function(row) {
      non_zero_values <- row[row != 0]
      if (length(non_zero_values) > 0) {
        return(mean(non_zero_values))
      } else {
        return(0)
      }
    })
  } else if (strict == T) {
    glm_coef <- apply(coef_df, 1, function(row) {
      if (any(row == 0)) {
        return(0)
      } else {
        return(mean(row))
      }
    })
  }
  
  glm_coef
}


find_glm_coef2 <- function(train_df,
                           pecked_part = 0.25,
                           n_of_pecking = 5,
                           lasso = FALSE,
                           strict = FALSE,
                           keep_original_1 = TRUE,
                           clustering_method = c("dbscan", "kmeans"),
                           dbscan_minPts = 5,
                           auto_eps_k = 4,
                           dbscan_noise = c("class0", "drop")) {
  
  clustering_method <- match.arg(clustering_method)
  dbscan_noise <- match.arg(dbscan_noise)
  
  for (k in 1:n_of_pecking) {
    set.seed(k)
    q <- pecked_part
    
    S1_ind <- which(train_df$S == 1)
    S0_ind <- which(train_df$S == 0)
    
    S1_ind_sample <- sample(S1_ind, length(S1_ind) * q)
    
    C1_minus_q <- train_df[setdiff(S1_ind, S1_ind_sample), ]
    C0_plus_q <- train_df[c(S0_ind, S1_ind_sample), ]
    
    C0_plus_q <- C0_plus_q %>%
      mutate_if(~ length(unique(.)) == 2, as.factor) %>%
      mutate(S = as.numeric(S) - 1)
    
    X <- select(C0_plus_q, -any_of(c("S", "Y")))
    X_scaled <- scale(X)
    
    
    # Klasteryzacja
    current_method <- clustering_method
    cluster_res <- NULL
    
    if (current_method == "dbscan") {
      kNNd <- dbscan::kNNdist(X_scaled, k = auto_eps_k)
      dists <- sort(kNNd)
      d2 <- diff(diff(dists))
      elbow <- which.max(d2)
      eps_final <- dists[elbow + 1]
      
      cluster_res <- dbscan::dbscan(X_scaled, eps = eps_final, minPts = dbscan_minPts)
      C0_plus_q$cluster <- cluster_res$cluster
      
      if (dbscan_noise == "drop") {
        C0_plus_q <- C0_plus_q %>% filter(cluster != 0)
      } else if (dbscan_noise == "class0") {
        C0_plus_q$cluster[C0_plus_q$cluster == 0] <- -1
      }
      
      if (length(unique(C0_plus_q$cluster)) < 2) {
        message(glue::glue("Iteration {k}: DBSCAN returned <2 clusters, switching to KMeans."))
        current_method <- "kmeans"
        cluster_res <- kmeans(X_scaled, centers = 2, nstart = 5)
        C0_plus_q$cluster <- as.numeric(cluster_res$cluster)
      }
    }
    
    if (current_method == "kmeans") {
      if (is.null(cluster_res)) {
        cluster_res <- kmeans(X_scaled, centers = 2, nstart = 5)
        C0_plus_q$cluster <- as.numeric(cluster_res$cluster)
      }
    }
    
    # Wybór klastra podobnego do S == 1
    like1clust <- C0_plus_q %>%
      group_by(cluster) %>%
      summarise(frac = mean(S)) %>%
      arrange(desc(frac)) %>%
      slice(1) %>%
      pull(cluster)
    
    C1 <- C0_plus_q %>% filter(cluster == like1clust) %>% select(-cluster)
    C0 <- C0_plus_q %>% filter(cluster != like1clust) %>% select(-cluster)
    
    if (dbscan_noise == "class0" && any(C0_plus_q$cluster == -1)) {
      noise_rows <- C0_plus_q %>% filter(cluster == -1) %>% select(-cluster)
      C0 <- bind_rows(C0, noise_rows)
    }
    
    if (!keep_original_1) {
      to_Y_modelling <- rbind.data.frame(
        cbind.data.frame(rbind.data.frame(C1_minus_q, C1), YY = 1),
        cbind.data.frame(C0, YY = 0)
      ) %>% select(-c(S)) %>% mutate_all(as.numeric)
    } else {
      to_Y_modelling <- rbind.data.frame(
        cbind.data.frame(rbind.data.frame(C1_minus_q, C1, C0 %>% filter(S == 1)), YY = 1),
        cbind.data.frame(C0 %>% filter(S == 0), YY = 0)
      ) %>% select(-c(S)) %>% mutate_all(as.numeric)
    }
    
    if (lasso) {
      x_to_lasso <- to_Y_modelling %>% select(-YY) %>% as.matrix()
      obj3 <- glmnet::cv.glmnet(x = x_to_lasso,
                                y = to_Y_modelling$YY,
                                standardize = TRUE,
                                intercept = TRUE,
                                family = "binomial")
      lambda <- obj3$lambda.min
      delta <- lambda * 0.5
      betasx1 <- coefficients(obj3, s = lambda)
      lasso_betas <- as.vector(betasx1)
      names(lasso_betas) <- row.names(betasx1)
      current_coef <- ifelse(abs(lasso_betas) > delta, lasso_betas, 0)
    } else {
      glm_clust <- glm(YY ~ ., data = to_Y_modelling, family = binomial)
      current_coef <- coef(glm_clust)
    }
    
    if (k == 1) {
      coef_df <- data.frame(row.names = names(current_coef))
      coef_df[[k]] <- current_coef
    } else {
      coef_df[[k]] <- current_coef
    }
  }
  
  if (!exists("coef_df")) {
    stop("No valid iterations completed (e.g., all DBSCAN runs failed and KMeans fallback not triggered).")
  }
  
  if (!strict) {
    glm_coef <- apply(coef_df, 1, function(row) {
      non_zero_values <- row[row != 0]
      if (length(non_zero_values) > 0) {
        return(mean(non_zero_values))
      } else {
        return(0)
      }
    })
  } else {
    glm_coef <- apply(coef_df, 1, function(row) {
      if (any(row == 0)) {
        return(0)
      } else {
        return(mean(row))
      }
    })
  }
  attr(glm_coef, "used_clustering_method") <- current_method
  return(glm_coef)
}

spy_sem <- function(df,
                    spy_frac = 0.2,
                    noise_level = 0.15,
                    em1_iter = 5,
                    em2_iter = 10,
                    verbose = TRUE) {
  
  library(dplyr)
  library(naivebayes)
  
  feature_cols <- setdiff(names(df), c("S","Y"))
  
  ### ===== Podział =====
  P <- df %>% filter(S == 1)
  M <- df %>% filter(S == 0)
  
  ### ===== Spies =====
  n_spy <- ceiling(spy_frac * nrow(P))
  spy_idx <- sample(nrow(P), n_spy)
  
  spies   <- P[spy_idx, ]
  P_clean <- P[-spy_idx, ]
  
  spies$spy    <- TRUE
  M$spy        <- FALSE
  P_clean$spy  <- FALSE
  
  MS <- bind_rows(M, spies)
  
  ### ====== FAZA 1: I-EM ======
  
  MS$Pr_pos      <- 0
  P_clean$Pr_pos <- 1
  
  for (i in 1:em1_iter) {
    
    train <- bind_rows(
      P_clean %>% mutate(Class = "pos", weight = 1),
      MS %>% mutate(Class = "neg", weight = 1 - Pr_pos)
    )
    
    model_em1 <- naive_bayes(
      x = train[, feature_cols],
      y = factor(train$Class),
      weights = train$weight
    )
    
    probs <- predict(model_em1, MS[, feature_cols], type = "prob")
    MS$Pr_pos <- probs[, "pos"]
    
    if (verbose) cat("EM1 iter:", i, "\n")
  }
  
  ### ====== Próg ze spies ======
  
  spy_probs <- MS %>% filter(spy) %>% pull(Pr_pos)
  t <- sort(spy_probs)[ceiling(noise_level * length(spy_probs))]
  
  N <- MS %>% filter(!spy & Pr_pos < t)
  U <- MS %>% filter(!spy & Pr_pos >= t)
  
  if (verbose) {
    cat("Likely negatives:", nrow(N), "\n")
    cat("Unlabeled:", nrow(U), "\n")
  }
  
  ### ====== FAZA 2: S-EM ======
  
  if (nrow(N) > 0) N$Pr_pos <- 0
  if (nrow(U) > 0) U$Pr_pos <- 0.5
  
  P_full <- bind_rows(P_clean, spies)
  P_full$Pr_pos <- 1
  
  all_data <- bind_rows(P_full, N, U)
  
  ### Flaga powodzenia EM2
  em2_success <- TRUE
  
  for (i in 1:em2_iter) {
    
    train <- all_data %>%
      mutate(
        Class  = ifelse(Pr_pos >= 0.5, "pos", "neg"),
        weight = ifelse(Class == "pos", Pr_pos, 1 - Pr_pos)
      ) %>%
      filter(weight > 0)
    
    # brak danych do trenowania
    if (nrow(train) == 0) {
      if (verbose) cat("EM2 iter:", i, "- empty train, stopping EM2\n")
      em2_success <- FALSE
      break
    }
    
    # tylko jedna klasa
    if (length(unique(train$Class)) < 2) {
      if (verbose) cat("EM2 iter:", i, "- skipping (one class)\n")
      em2_success <- FALSE
      break
    }
    
    model_em2 <- naive_bayes(
      x = train[, feature_cols],
      y = factor(train$Class),
      weights = train$weight
    )
    
    probs <- predict(model_em2, all_data[, feature_cols], type = "prob")
    all_data$Pr_pos <- probs[, "pos"]
    
    if (verbose) cat("EM2 iter:", i, "\n")
  }
  
  ### ===== Fallback do EM1 =====
  
  if (!em2_success) {
    if (verbose) cat("EM2 failed — returning EM1 model\n")
    
    all_data <- bind_rows(P_clean, spies, N, U)
    all_data$final_class <- ifelse(all_data$Pr_pos > 0.5, 1, 0)
    
    return(list(
      model      = model_em1,
      model_type = "em1",
      data       = all_data,
      threshold  = t,
      N          = N,
      U          = U
    ))
  }
  
  ### ===== Wynik z EM2 =====
  
  all_data$final_class <- ifelse(all_data$Pr_pos > 0.5, 1, 0)
  
  list(
    model      = model_em2,
    model_type = "em2",
    data       = all_data,
    threshold  = t,
    N          = N,
    U          = U
  )
}

round_to_set <- function(x, set = c(0.3, 0.5, 0.8)) {
  set[apply(abs(outer(x, set, "-")), 1, which.min)]
}
