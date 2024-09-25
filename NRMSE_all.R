#NRMSE for original dataset (Groupwise data normalization )
nrmse_eval <- function (data_input, groups){
  suppressMessages({
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
  
  #VSN Normalization function
  VSN_Norm <- function(dat) {
    dat<-as.data.frame(dat)
    vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
    colnames(vsnNormed) <- colnames(dat)
    row.names(vsnNormed) <- rownames(dat)
    return(as.matrix(vsnNormed))
  }
  
  #VSN normalized data
  vsn.dat <- do.call("cbind", lapply(group_data, VSN_Norm))
  
  #Loess normalization
  LOESS_Norm <- function(dat) {
    newdata<- as.data.frame(log2(dat))
    cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(newdata), method="fast")
    colnames(cycLoessNormed) <- colnames(newdata)
    row.names(cycLoessNormed) <- rownames(newdata)
    return(as.matrix(cycLoessNormed))
  }
  
  
  #Loess normalized data
  loess.dat <- do.call("cbind", lapply(group_data, LOESS_Norm))
  
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
  
  #RLR normalized data
  rlr.dat <- do.call("cbind", lapply(group_data, RLR_Norm))
  
  #Grouping of normalized datasets
  vsn_group_data <- grouping_data(vsn.dat, groups)
  
  loess_group_data <- grouping_data(loess.dat, groups)
  
  rlr_group_data <- grouping_data(rlr.dat, groups)
  
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
  
  #VSN and KNN
  vsn.knn.dat <- do.call("cbind", lapply(vsn_group_data, KNN_Imputation))
  
  #VSN and LLS
  vsn.lls.dat <- do.call("cbind", lapply(vsn_group_data, LLS_Imputation))
  
  #VSN and SVD
  vsn.svd.dat <- do.call("cbind", lapply(vsn_group_data, SVD_Imputation))
  
  #LOESS and KNN
  loess.knn.dat <- do.call("cbind", lapply(loess_group_data, KNN_Imputation))
  
  #LOESS and LLS
  loess.lls.dat <- do.call("cbind", lapply(loess_group_data, LLS_Imputation))
  
  #LOESS and SVD
  loess.svd.dat <- do.call("cbind", lapply(loess_group_data, SVD_Imputation))
  
  #RLR and KNN
  rlr.knn.dat <- do.call("cbind", lapply(rlr_group_data, KNN_Imputation))
  
  #RLR and LLS
  rlr.lls.dat <- do.call("cbind", lapply(rlr_group_data, LLS_Imputation))
  
  #RLR and SVD
  rlr.svd.dat <- do.call("cbind", lapply(rlr_group_data, SVD_Imputation))
  
  
  # Apply imputation methods and handle missing values
  vsn.knn.dat[is.na(vsn.knn.dat)] <- loess.knn.dat[is.na(loess.knn.dat)] <- rlr.knn.dat[is.na(rlr.knn.dat)] <- 0
  vsn.lls.dat[is.na(vsn.lls.dat)] <- loess.lls.dat[is.na(loess.lls.dat)] <- rlr.lls.dat[is.na(rlr.lls.dat)] <- 0
  vsn.svd.dat[is.na(vsn.svd.dat)] <- loess.svd.dat[is.na(loess.svd.dat)] <- rlr.svd.dat[is.na(rlr.svd.dat)] <- 0
  
  #Transposing the sample to wide format
  new_sample <- as.data.frame(t(groups))
  names(new_sample) <- new_sample[1,]
  sample <- new_sample[-1,]
  
  #Changing all data files column names
  new_colnames <- colnames(sample)
  colnames(vsn.knn.dat) <- new_colnames
  colnames(vsn.lls.dat) <- new_colnames
  colnames(vsn.svd.dat) <- new_colnames
  colnames(loess.knn.dat) <- new_colnames
  colnames(loess.lls.dat) <- new_colnames
  colnames(loess.svd.dat) <- new_colnames
  colnames(rlr.knn.dat) <- new_colnames
  colnames(rlr.lls.dat) <- new_colnames
  colnames(rlr.svd.dat) <- new_colnames
  
  com_data2[is.na(com_data2)] <- 0

  nrmse <- function(ximp, xtrue) {
    # Convert both inputs to numeric vectors
    ximp_vector <- as.numeric(as.matrix(ximp))
    xtrue_vector <- as.numeric(as.matrix(xtrue))
    
    # Calculate NRMSE
    sqrt(mean((ximp_vector - xtrue_vector)^2) / stats::var(xtrue_vector))
  }
  
  nrmse_vsn_knn <- nrmse(2^vsn.knn.dat, com_data2)
  nrmse_vsn_lls <- nrmse(2^vsn.lls.dat, com_data2)
  nrmse_vsn_svd <- nrmse(2^vsn.svd.dat, com_data2)
  nrmse_loess_knn <- nrmse(2^loess.knn.dat, com_data2)
  nrmse_loess_lls <- nrmse(2^loess.lls.dat, com_data2)
  nrmse_loess_svd <- nrmse(2^loess.svd.dat, com_data2)
  nrmse_rlr_knn <- nrmse(2^rlr.knn.dat, com_data2)
  nrmse_rlr_lls <- nrmse(2^rlr.lls.dat, com_data2)
  nrmse_rlr_svd <- nrmse(2^rlr.svd.dat, com_data2)
  
  nrmse <- c(nrmse_vsn_knn, nrmse_vsn_lls, nrmse_vsn_svd,
             nrmse_loess_knn, nrmse_loess_lls, nrmse_loess_svd,
             nrmse_rlr_knn, nrmse_rlr_lls, nrmse_rlr_svd)
  
  nrmse_1 <- as.data.frame(round(nrmse, digits = 5))
  
  nrmse_result <- cbind(Combination = c("vsn_knn", "vsn_lls", "vsn_svd",
                                        "loess_knn", "loess_lls", "loess_svd",
                                        "rlr_knn", "rlr_lls", "rlr_svd"), nrmse_1)
  colnames(nrmse_result) <- c("Combination", "NRMSE")  
  
  result <- list("NRMSE result" = nrmse_result)
  
  return(result)
  })
})
}

#NRMSE for original dataset (Complete data normalization )
nrmse_comp <- function (data_input, groups){
  suppressMessages({
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
  
  vsn.dat <- VSN_Norm(data_input)
  
  # Loess normalization function
  LOESS_Norm <- function(dat) {
    newdata <- as.data.frame(log2(dat))
    cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(newdata), method = "fast")
    colnames(cycLoessNormed) <- colnames(newdata)
    row.names(cycLoessNormed) <- rownames(newdata)
    return(as.matrix(cycLoessNormed))
  }
  
  loess.dat <- LOESS_Norm(data_input)
  
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
  
  rlr.dat <- RLR_Norm(data_input)
  
  
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
  
  # Apply imputation methods
  vsn.knn.dat <- KNN_Imputation(vsn.dat)
  vsn.lls.dat <- LLS_Imputation(vsn.dat)
  vsn.svd.dat <- SVD_Imputation(vsn.dat)
  
  loess.knn.dat <- KNN_Imputation(loess.dat)
  loess.lls.dat <- LLS_Imputation(loess.dat)
  loess.svd.dat <- SVD_Imputation(loess.dat)
  
  rlr.knn.dat <- KNN_Imputation(rlr.dat)
  rlr.lls.dat <- LLS_Imputation(rlr.dat)
  rlr.svd.dat <- SVD_Imputation(rlr.dat)
  
  # Apply imputation methods and handle missing values
  vsn.knn.dat[is.na(vsn.knn.dat)] <- loess.knn.dat[is.na(rlr.knn.dat)] <- rlr.knn.dat[is.na(rlr.knn.dat)] <- 0
  vsn.lls.dat[is.na(vsn.lls.dat)] <- loess.lls.dat[is.na(rlr.lls.dat)] <- rlr.lls.dat[is.na(rlr.lls.dat)] <- 0
  vsn.svd.dat[is.na(vsn.svd.dat)] <- loess.svd.dat[is.na(rlr.svd.dat)] <- rlr.svd.dat[is.na(rlr.svd.dat)] <- 0
  
  nrmse <- function(ximp, xtrue) {
    # Convert both inputs to numeric vectors
    ximp_vector <- as.numeric(as.matrix(ximp))
    xtrue_vector <- as.numeric(as.matrix(xtrue))
    
    # Calculate NRMSE
    sqrt(mean((ximp_vector - xtrue_vector)^2) / stats::var(xtrue_vector))
  }
  
  data_input[is.na(data_input)] <- 0
  
  nrmse_vsn_knn <- nrmse(2^vsn.knn.dat, data_input)
  nrmse_vsn_lls <- nrmse(2^vsn.lls.dat, data_input)
  nrmse_vsn_svd <- nrmse(2^vsn.svd.dat, data_input)
  nrmse_loess_knn <- nrmse(2^loess.knn.dat, data_input)
  nrmse_loess_lls <- nrmse(2^loess.lls.dat, data_input)
  nrmse_loess_svd <- nrmse(2^loess.svd.dat, data_input)
  nrmse_rlr_knn <- nrmse(2^rlr.knn.dat, data_input)
  nrmse_rlr_lls <- nrmse(2^rlr.lls.dat, data_input)
  nrmse_rlr_svd <- nrmse(2^rlr.svd.dat, data_input)
  
  nrmse <- c(nrmse_vsn_knn, nrmse_vsn_lls, nrmse_vsn_svd,
             nrmse_loess_knn, nrmse_loess_lls, nrmse_loess_svd,
             nrmse_rlr_knn, nrmse_rlr_lls, nrmse_rlr_svd)
  
  nrmse_1 <- as.data.frame(round(nrmse, digits = 5))
  
  nrmse_result <- cbind(Combination = c("vsn_knn", "vsn_lls", "vsn_svd",
                                        "loess_knn", "loess_lls", "loess_svd",
                                        "rlr_knn", "rlr_lls", "rlr_svd"), nrmse_1)
  colnames(nrmse_result) <- c("Combination", "NRMSE")  
  
  result <- list("NRMSE result" = nrmse_result)
  
  return(result)
    })
  })
}


