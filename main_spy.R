source('_helpers.R')
library(pROC)
library(dplyr)
library(caret)
library(pROC)
library(MLmetrics)

datasets <- load_datasets(more_then_n_unique_values = 5, dataset_type = 'n', majority_class = 1)
list2env(datasets, envir = .GlobalEnv)


#if you have dowloaded datasets you can uncomment below
# files <- list.files("datasets", full.names = TRUE, pattern = "\\.rds$")
# 
# for (f in files) {
#   name <- tools::file_path_sans_ext(basename(f))
#   assign(name, readRDS(f))
# }

dsets<-
  ls()[sapply(ls(), function(x) class(get(x))) == 'data.frame']

for (i in dsets) {
  
  current_frame <- get(i)

  experi_max1 <- ifelse(ncol(current_frame)<100,100,10)
  experi_max2 <- ifelse(nrow(current_frame)<10000,100,10)
  
  experi_max <- min(experi_max1,experi_max2)
  
  for (num in 1:experi_max) {
    set.seed(num)
    gc()
    for (c_calc in c(0.3,0.5,0.8)) {
      
      #current_frame_with_S <- non_scar_labelling.mvc(current_frame,c_calc, n_vars = 2)
      current_frame_with_S <- scar_labelling(current_frame,c_calc)
      index <- createDataPartition(current_frame_with_S$S, p = 0.7, list = F)
      train <- current_frame_with_S[index,]
      test <- current_frame_with_S[-index,]
      
      feature_cols <- setdiff(names(train), c("S","Y"))
      
      {res <- spy_sem(train,
                     spy_frac = 0.2,
                     noise_level = 0.15,
                     em1_iter = 5,
                     em2_iter = 10)} %>% system.time(.) -> time
      spy_time <- time["elapsed"] %>% unname
      
      probs <- predict(res$model, test[, feature_cols], type = "prob")
      
      Y_pred <- probs[,"pos"]
      Y_pred_label <- ifelse(probs[,"pos"]>0.5,1,0)
      y_test <- test$Y
      
      c_calculated <- sum(current_frame_with_S$S)/sum(current_frame_with_S$Y)
      
      res_row <- data.frame(
        ts = Sys.time(),
        df = i,
        experi_num = num,
        calc_c = round(c_calculated,2),
        spy_auc = suppressMessages(round(auc(y_test, Y_pred),3)),
        spy_acc = suppressMessages(round(Accuracy(Y_pred_label, y_test),3)),
        spy_f1 = suppressMessages(round(F1_Score(Y_pred_label, y_test),3)),
        spy_time = spy_time
        
        )
      
      if (!exists('res_df')) {
        res_df <- res_row
      } else {
        res_df <- rbind.data.frame(res_df,res_row)
      }
      
      print(res_row)
      
      
    }
  }
  
  short_commit_hash <- system("git rev-parse --short HEAD", intern = TRUE)
  output_filename <- paste0('res_scar_spy_',short_commit_hash,'.RDS')
  saveRDS(object = res_df, file = output_filename)
}

