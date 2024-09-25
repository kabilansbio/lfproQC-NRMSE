#NRMSE imputation followed by normalization in original dataset (Groupwise data normalization )
nrmse_eval_imp_norm <- function (data_input, groups){
  suppressWarnings({
  
  # Converting all zeros to NAs
  data_input[data_input == 0] <- NA
  
  complete_data_fn <- function(data, groups) {
    # Rename the columns in the groups data
    colnames(groups) <- c("name", "group")
    
    # Add a new column to preserve the original order of rows
    data <- data %>% dplyr::mutate(original_order = dplyr::row_number())
    
    # Rename the first column of data
    colnames(data)[1] <- "rowid"
    
    # Reshape the data into long format
    long_data <- data %>%
      tidyr::pivot_longer(-c(rowid, original_order), names_to = "name", values_to = "mass")
    
    # Merge with groups to add group information
    long_data <- long_data %>%
      dplyr::left_join(groups, by = "name")
    
    # Identify rows where any group has all missing values
    group_summary <- long_data %>%
      dplyr::group_by(rowid, group) %>%
      dplyr::summarise(all_na = all(is.na(mass)), .groups = 'drop') %>%
      dplyr::ungroup()
    
    # Identify rows to remove (where any group has all missing values)
    rows_to_remove <- group_summary %>%
      dplyr::group_by(rowid) %>%
      dplyr::summarise(remove = any(all_na)) %>%
      dplyr::filter(remove) %>%
      dplyr::pull(rowid) %>%
      unique()
    
    # Filter out rows with any completely missing group
    filtered_data <- long_data %>%
      dplyr::filter(!(rowid %in% rows_to_remove)) %>%
      dplyr::select(-group)
    
    # Reshape back to wide format and reorder based on the original order
    com_data <- filtered_data %>%
      tidyr::pivot_wider(names_from = name, values_from = mass) %>%
      dplyr::arrange(original_order) %>%
      dplyr::select(-original_order)
    
    return(com_data)
  } 
  
  com_data <- complete_data_fn (data_input, groups)
  
  #Giving original name to the first column
  com_data1 <- com_data
  colnames(com_data1)[1] <- colnames(data_input)[1]
  
  #Extracting the first column contains ID information
  com_data_ID <- com_data1[,1]
  
  #Removing the ID column and selecting remaining data
  com_data2 <- com_data1[,-1]
  
  #Giving original name to the first column
  com_data1 <- com_data
  colnames(com_data1)[1] <- colnames(data_input)[1]
  
  #Extracting the first column contains ID information
  com_data_ID <- com_data1[,1]
  
  #Removing the ID column and selecting remaining data
  com_data2 <- com_data1[,-1]
  #Grouping of dataframe according to sample groups
  grouping_data <- function(df, groups) {  # df = dataframe, groups = groups dataframe
    # Rename columns in groups dataframe for consistency
    colnames(groups) <- c("name", "group")
    
    # Create a list to store grouped columns
    grouped_data <- list()
    
    # Loop through each unique group
    for (group_name in unique(groups$group)) {
      # Get the names of the columns that belong to the current group
      group_columns <- groups$name[groups$group == group_name]
      
      # Extract the columns from df that match the current group columns
      grouped_data[[group_name]] <- df[, group_columns, drop = FALSE]
    }
    
    return(grouped_data)
  }
  
  #Grouping of dataframe as a triplicate groups
  group_data <- grouping_data(com_data2, groups)
  
  #Imputation of normalized datasets
  #KNN imputation
  KNN_Imputation <- function (dat)
  {
    resultkNN <- VIM::kNN(dat, numFun = laeken::weightedMean, weightDist = TRUE,
                          imp_var = FALSE, k= 10)
    return(resultkNN)
  }
  
  #LLS imputation
  LLS_Imputation <- function (dat)
  {
    resultLLS <- pcaMethods::llsImpute(dat, k=2, correlation = "pearson", allVariables = TRUE)
    dataSet.imputed <- resultLLS@completeObs
    return(dataSet.imputed)
  }
  
  #SVD imputation
  SVD_Imputation <- function (dat)
  {
    resultSVD <- pcaMethods::pca(dat, method = "svdImpute", nPcs = 2)
    dataSet.imputed <- resultSVD@completeObs
    return(dataSet.imputed)
  }
  
  #KNN 
  knn.dat <- do.call("cbind", lapply(group_data, KNN_Imputation))
  colnames(knn.dat) <- colnames(com_data2)
  knn_group_data <- grouping_data(knn.dat, groups)
  
  #LLS 
  lls.dat <- do.call("cbind", lapply(group_data, LLS_Imputation))
  colnames(lls.dat) <- colnames(com_data2)
  lls_group_data <- grouping_data(lls.dat, groups)
  
  #SVD 
  svd.dat <- do.call("cbind", lapply(group_data, SVD_Imputation))
  colnames(svd.dat) <- colnames(com_data2)
  svd_group_data <- grouping_data(svd.dat, groups)
  
  #VSN Normalization function
  VSN_Norm <- function(dat) {
    dat<-as.data.frame(dat)
    vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
    colnames(vsnNormed) <- colnames(dat)
    row.names(vsnNormed) <- rownames(dat)
    return(as.matrix(vsnNormed))
  }
  
  #Loess normalization
  LOESS_Norm <- function(dat) {
    newdata<- as.data.frame(log2(dat))
    cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(newdata), method="fast")
    colnames(cycLoessNormed) <- colnames(newdata)
    row.names(cycLoessNormed) <- rownames(newdata)
    return(as.matrix(cycLoessNormed))
  }
  
  
  #RLR normalization
  RLR_Norm <- function (dat) {
    log2Matrix <- log2(dat)
    log2Matrix <- as.matrix(log2Matrix)
    
    # Extract the numeric component from the list
    log2Matrix <- unlist(log2Matrix)
    log2Matrix[is.infinite(log2Matrix)] <- 0
    
    sampleLog2Median <- matrixStats::rowMedians(log2Matrix, na.rm = TRUE)
    calculateRLMForCol <- function(colIndex, sampleLog2Median,
                                   log2Matrix) {
      lrFit <- MASS::rlm(as.matrix(log2Matrix[, colIndex]) ~
                           sampleLog2Median, na.action = stats::na.exclude, maxit = 20)
      coeffs <- lrFit$coefficients
      coefIntercept <- coeffs[1]
      coefSlope <- coeffs[2]
      globalFittedRLRCol <- (log2Matrix[, colIndex] - coefIntercept)/coefSlope
      globalFittedRLRCol
    }
    
    data <- log2Matrix
    globalFittedRLR <- vapply(seq_len(ncol(data)), calculateRLMForCol,
                              rep(0, nrow(data)), sampleLog2Median = sampleLog2Median,
                              log2Matrix = data)
    colnames(globalFittedRLR) <- colnames(dat)
    globalFittedRLR[globalFittedRLR < 0] <- 0
    globalFittedRLR
  }

  #KNN and VSN
  knn.vsn.dat <- do.call("cbind", lapply(knn_group_data, VSN_Norm))
  
  #KNN and LOESS
  knn.loess.dat <- do.call("cbind", lapply(knn_group_data, LOESS_Norm))
  
  #KNN and RLR
  knn.rlr.dat <- do.call("cbind", lapply(knn_group_data, RLR_Norm))
  
  #LLS and VSN
  lls.vsn.dat <- do.call("cbind", lapply(lls_group_data, VSN_Norm))
  
  #LLS and LOESS
  lls.loess.dat <- do.call("cbind", lapply(lls_group_data, LOESS_Norm))
  
  #LLS and RLR
  lls.rlr.dat <- do.call("cbind", lapply(lls_group_data, RLR_Norm))
  
  #SVD and VSN
  svd.vsn.dat <- do.call("cbind", lapply(svd_group_data, VSN_Norm))
  
  #SVD and LOESS
  svd.loess.dat <- do.call("cbind", lapply(svd_group_data, LOESS_Norm))
  
  #SVD and RLR
  svd.rlr.dat <- do.call("cbind", lapply(svd_group_data, RLR_Norm))

  knn.vsn.dat[is.na(knn.vsn.dat)] <- knn.loess.dat[is.na(knn.loess.dat)] <- knn.rlr.dat[is.na(knn.rlr.dat)] <-0
  lls.vsn.dat[is.na(lls.vsn.dat)] <- lls.loess.dat[is.na(lls.loess.dat)] <- lls.rlr.dat[is.na(lls.rlr.dat)] <-0
  svd.vsn.dat[is.na(svd.vsn.dat)] <- svd.loess.dat[is.na(svd.loess.dat)] <- svd.rlr.dat[is.na(svd.rlr.dat)] <-0
  
  
  #Transposing the sample to wide format
  new_sample <- as.data.frame(t(groups))
  names(new_sample) <- new_sample[1,]
  sample <- new_sample[-1,]
  
  #Changing all data files column names
  new_colnames <- colnames(sample)
  colnames(knn.vsn.dat) <- new_colnames
  colnames(knn.loess.dat) <- new_colnames
  colnames(knn.rlr.dat) <- new_colnames
  colnames(lls.vsn.dat) <- new_colnames
  colnames(lls.loess.dat) <- new_colnames
  colnames(lls.rlr.dat) <- new_colnames
  colnames(svd.vsn.dat) <- new_colnames
  colnames(svd.loess.dat) <- new_colnames
  colnames(svd.rlr.dat) <- new_colnames
  
  com_data2[is.na(com_data2)] <- 0
  
  nrmse <- function(ximp, xtrue) {
    # Convert both inputs to numeric vectors
    ximp_vector <- as.numeric(as.matrix(ximp))
    xtrue_vector <- as.numeric(as.matrix(xtrue))
    
    # Ensure the vectors have the same length
    if (length(ximp_vector) != length(xtrue_vector)) {
      stop("The inputs must have the same length.")
    }
    
    # Calculate NRMSE
    sqrt(mean((ximp_vector - xtrue_vector)^2) / stats::var(xtrue_vector))
  }
  
  # Calculate normalized root mean squared error (NRMSE)
  nrmse_knn_vsn <- nrmse(2^knn.vsn.dat, com_data2)
  nrmse_knn_loess <- nrmse(2^knn.loess.dat, com_data2)
  nrmse_knn_rlr <- nrmse(2^knn.rlr.dat, com_data2)
  nrmse_lls_vsn <- nrmse(2^lls.vsn.dat, com_data2)
  nrmse_lls_loess <- nrmse(2^lls.loess.dat, com_data2)
  nrmse_lls_rlr <- nrmse(2^lls.rlr.dat, com_data2)
  nrmse_svd_vsn <- nrmse(2^svd.vsn.dat, com_data2)
  nrmse_svd_loess <- nrmse(2^svd.loess.dat, com_data2)
  nrmse_svd_rlr <- nrmse(2^svd.rlr.dat, com_data2)
  
  
  # Create a vector with the updated variable names
  nrmse <- c(nrmse_knn_vsn, nrmse_knn_loess, nrmse_knn_rlr,
             nrmse_lls_vsn, nrmse_lls_loess, nrmse_lls_rlr,
             nrmse_svd_vsn, nrmse_svd_loess, nrmse_svd_rlr)
  
  nrmse_1 <- as.data.frame(round(nrmse, digits = 5))
  
  nrmse_result <- cbind(Combination =  c("knn_vsn", "knn_loess", "knn_rlr",
                                         "lls_vsn", "lls_loess", "lls_rlr",
                                         "svd_vsn", "svd_loess", "svd_rlr"), nrmse_1)
  
  colnames(nrmse_result) <- c("Combination", "NRMSE")  
  
  result <- list("NRMSE result" = nrmse_result)
  
  return(result)
})
}

#NRMSE imputation followed by normalization in original dataset (Complete data normalization )
nrmse_comp_imp_norm  <- function (data_input, groups){
  suppressWarnings({
  
  # Converting all zeros to NAs
  data_input[data_input == 0] <- NA
  
  complete_data_fn <- function(data, groups) {
    # Rename the columns in the groups data
    colnames(groups) <- c("name", "group")
    
    # Add a new column to preserve the original order of rows
    data <- data %>% dplyr::mutate(original_order = dplyr::row_number())
    
    # Rename the first column of data
    colnames(data)[1] <- "rowid"
    
    # Reshape the data into long format
    long_data <- data %>%
      tidyr::pivot_longer(-c(rowid, original_order), names_to = "name", values_to = "mass")
    
    # Merge with groups to add group information
    long_data <- long_data %>%
      dplyr::left_join(groups, by = "name")
    
    # Identify rows where any group has all missing values
    group_summary <- long_data %>%
      dplyr::group_by(rowid, group) %>%
      dplyr::summarise(all_na = all(is.na(mass)), .groups = 'drop') %>%
      dplyr::ungroup()
    
    # Identify rows to remove (where any group has all missing values)
    rows_to_remove <- group_summary %>%
      dplyr::group_by(rowid) %>%
      dplyr::summarise(remove = any(all_na)) %>%
      dplyr::filter(remove) %>%
      dplyr::pull(rowid) %>%
      unique()
    
    # Filter out rows with any completely missing group
    filtered_data <- long_data %>%
      dplyr::filter(!(rowid %in% rows_to_remove)) %>%
      dplyr::select(-group)
    
    # Reshape back to wide format and reorder based on the original order
    com_data <- filtered_data %>%
      tidyr::pivot_wider(names_from = name, values_from = mass) %>%
      dplyr::arrange(original_order) %>%
      dplyr::select(-original_order)
    
    return(com_data)
  } 
  
  com_data <- complete_data_fn (data_input, groups)
  
  #Giving original name to the first column
  com_data1 <- com_data
  colnames(com_data1)[1] <- colnames(data_input)[1]
  
  #Extracting the first column contains ID information
  com_data_ID <- com_data1[,1]
  
  #Removing the ID column and selecting remaining data
  com_data2 <- com_data1[,-1]
  
  #Giving original name to the first column
  com_data1 <- com_data
  colnames(com_data1)[1] <- colnames(data_input)[1]
  
  #Extracting the first column contains ID information
  com_data_ID <- com_data1[,1]
  
  #Removing the ID column and selecting remaining data
  com_data2 <- com_data1[,-1]
  
  data_input <- as.data.frame(com_data2)
  
  # VSN Normalization function
  VSN_Norm <- function(dat) {
    dat <- as.data.frame(dat)
    vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
    colnames(vsnNormed) <- colnames(dat)
    row.names(vsnNormed) <- rownames(dat)
    return(as.matrix(vsnNormed))
  }
  
  # Loess normalization function
  LOESS_Norm <- function(dat) {
    newdata <- as.data.frame(log2(dat))
    cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(newdata), method = "fast")
    colnames(cycLoessNormed) <- colnames(newdata)
    row.names(cycLoessNormed) <- rownames(newdata)
    return(as.matrix(cycLoessNormed))
  }
  
  # RLR normalization function
  RLR_Norm <- function(dat) {
    log2Matrix <- log2(dat)
    log2Matrix <- as.matrix(log2Matrix)
    log2Matrix[is.infinite(log2Matrix)] <- 0
    sampleLog2Median <- matrixStats::rowMedians(log2Matrix, na.rm = TRUE)
    
    calculateRLMForCol <- function(colIndex, sampleLog2Median, log2Matrix) {
      lrFit <- MASS::rlm(as.matrix(log2Matrix[, colIndex]) ~ sampleLog2Median, na.action = stats::na.exclude, maxit = 20)
      coeffs <- lrFit$coefficients
      coefIntercept <- coeffs[1]
      coefSlope <- coeffs[2]
      globalFittedRLRCol <- (log2Matrix[, colIndex] - coefIntercept) / coefSlope
      return(globalFittedRLRCol)
    }
    
    data <- log2Matrix
    globalFittedRLR <- vapply(seq_len(ncol(data)), calculateRLMForCol, rep(0, nrow(data)), sampleLog2Median = sampleLog2Median, log2Matrix = data)
    colnames(globalFittedRLR) <- colnames(dat)
    globalFittedRLR[globalFittedRLR < 0] <- 0
    return(globalFittedRLR)
  }
  
  KNN_Imputation <- function(dat) {
    resultkNN <- VIM::kNN(dat, numFun = laeken::weightedMean, weightDist = TRUE, imp_var = FALSE, k = 10)
    return(resultkNN)
  }
  
  LLS_Imputation <- function(dat) {
    resultLLS <- pcaMethods::llsImpute(dat, k = 2, correlation = "pearson", allVariables = TRUE)
    dataSet.imputed <- resultLLS@completeObs
    return(dataSet.imputed)
  }
  
  SVD_Imputation <- function(dat) {
    resultSVD <- pcaMethods::pca(dat, method = "svdImpute", nPcs = 2)
    dataSet.imputed <- resultSVD@completeObs
    return(dataSet.imputed)
  }
  
  #KNN and VSN
  knn.dat <- KNN_Imputation(data_input)
  
  #LLS and VSN
  lls.dat <- LLS_Imputation(data_input)
  
  #SVD and VSN
  svd.dat <- SVD_Imputation(data_input)
  
  # Apply normalization methods
  knn.vsn.dat <- VSN_Norm(knn.dat)
  knn.loess.dat <- LOESS_Norm(knn.dat)
  knn.rlr.dat <- RLR_Norm(knn.dat)
  
  lls.vsn.dat <- VSN_Norm(lls.dat)
  lls.loess.dat <- LOESS_Norm(lls.dat)
  lls.rlr.dat <- RLR_Norm(lls.dat)
  
  svd.vsn.dat <- VSN_Norm(svd.dat)
  svd.loess.dat <- LOESS_Norm(svd.dat)
  svd.rlr.dat <- RLR_Norm(svd.dat)
  
  knn.vsn.dat[is.na(knn.vsn.dat)] <- knn.loess.dat[is.na(knn.loess.dat)] <- knn.rlr.dat[is.na(knn.rlr.dat)] <-0
  lls.vsn.dat[is.na(lls.vsn.dat)] <- lls.loess.dat[is.na(lls.loess.dat)] <- lls.rlr.dat[is.na(lls.rlr.dat)] <-0
  svd.vsn.dat[is.na(svd.vsn.dat)] <- svd.loess.dat[is.na(svd.loess.dat)] <- svd.rlr.dat[is.na(svd.rlr.dat)] <-0
  
  nrmse <- function(ximp, xtrue) {
    # Convert both inputs to numeric vectors
    ximp_vector <- as.numeric(as.matrix(ximp))
    xtrue_vector <- as.numeric(as.matrix(xtrue))
    
    # Ensure the vectors have the same length
    if (length(ximp_vector) != length(xtrue_vector)) {
      stop("The inputs must have the same length.")
    }
    
    # Calculate NRMSE
    sqrt(mean((ximp_vector - xtrue_vector)^2) / stats::var(xtrue_vector))
  }
  
  data_input[is.na(data_input)] <- 0
  
  nrmse_knn_vsn <- nrmse(2^knn.vsn.dat, data_input)
  nrmse_knn_loess <- nrmse(2^knn.loess.dat, data_input)
  nrmse_knn_rlr <- nrmse(2^knn.rlr.dat, data_input)
  nrmse_lls_vsn <- nrmse(2^lls.vsn.dat, data_input)
  nrmse_lls_loess <- nrmse(2^lls.loess.dat, data_input)
  nrmse_lls_rlr <- nrmse(2^lls.rlr.dat, data_input)
  nrmse_svd_vsn <- nrmse(2^svd.vsn.dat, data_input)
  nrmse_svd_loess <- nrmse(2^svd.loess.dat, data_input)
  nrmse_svd_rlr <- nrmse(2^svd.rlr.dat, data_input)
  
  nrmse <- c(nrmse_knn_vsn, nrmse_knn_loess, nrmse_knn_rlr,
             nrmse_lls_vsn, nrmse_lls_loess, nrmse_lls_rlr,
             nrmse_svd_vsn, nrmse_svd_loess, nrmse_svd_rlr)
  
  nrmse_1 <- as.data.frame(round(nrmse, digits = 5))
  
  nrmse_result <- cbind(Combination = c("knn_vsn", "knn_loess", "knn_rlr",
                                        "lls_vsn", "lls_loess", "lls_rlr",
                                        "svd_vsn", "svd_loess", "svd_rlr"), nrmse_1)
  colnames(nrmse_result) <- c("Combination", "NRMSE")  
  
  result <- list("NRMSE result" = nrmse_result)
  
  return(result)
  })
}


