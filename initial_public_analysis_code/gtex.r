# The BiT-age method, one of the first methods we'll try, was designed and
# published by David Meyer and Bjorn Schumacher
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7963339/)
# (https://github.com/Meyer-DH/AgingClock)

the_glmnet_func <- function(df) {
  library(caret) # install.packages(c("caret", "glmnet"))
  library(glmnet)
  
  ages <- df$age
  df <- df[,-1]
  variances <- apply(df, 2, var)
  df <- df[, order(variances, decreasing = TRUE)[1:2000]]
  df$age <- ages
  
  # adjuvant subfunction
  get_best_result = function(caret_fit) { # from https://daviddalpiaz.github.io/r4sl/elastic-net.html
    best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
    best_result = caret_fit$results[best, ]
    rownames(best_result) = NULL
    best_result
  }
  
  set.seed(42)
  cv_3 = trainControl(method = "cv", number = 3)
  
  hit_elnet = train(
    age ~ ., data = df,
    method = "glmnet",
    trControl = cv_3
  )
  
  return(get_best_result(hit_elnet))
}

gtex_data <- readRDS("GTEX_21blood_samples.rds")
dim(gtex_data) # 21 54272
gtex_data <- gtex_data[,colSums(gtex_data)!=0] # remove genes that are zero-only
dim(gtex_data) # 21 39702

gtex_ages <- gtex_data[,"age"]
gtex_tpms <- gtex_data[,colnames(gtex_data) != "age"]
dim(gtex_tpms) # 21 39701

mod_tpms <- gtex_data
mod_tpms[mod_tpms==0] <- NA # This is justified by this: https://github.com/Meyer-DH/AgingClock/blob/b6980add65779be9da8af883409bed5f6f5e0e87/src/biological_age_prediction.py#L17

# Binarize based on whether TPM (similar to CPM) is above median or not.
sample_median_TPMs <- apply(mod_tpms, 1, FUN = function(x) {return(median(x, na.rm=TRUE))})
for (i in 1:nrow(mod_tpms)) {
  mod_tpms[i,] <- ifelse(mod_tpms[i,] > sample_median_TPMs[i], 1, 0)
}
mod_tpms[is.na(mod_tpms)] <- 0


# Doing the model with no binarisation ------
df <- as.data.frame(cbind(age=gtex_ages, gtex_tpms))
nobLM <- the_glmnet_func(df)
nobLM

# Doing the model with    binarisation ------
df <- as.data.frame(cbind(age=gtex_ages, mod_tpms))
binLM <- the_glmnet_func(df)
binLM

# Too few samples to be sure, but so far it looks like BiT-age is not a better predictor in blood samples from GTEX
