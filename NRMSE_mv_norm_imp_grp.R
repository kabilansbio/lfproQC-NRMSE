nrmse_group_mv_norm_imp <- function (data_input, groups){
  suppressMessages({
    suppressWarnings({
      
      # Store original row names
      original_row_names <- rownames(data_input)
      
      # Remove the first column while preserving the original row names
      data_input <- as.data.frame(data_input[,-1])
      rownames(data_input) <- original_row_names
      
      # Get the complete data (omit rows with any NA values)
      new_data <- stats::na.omit(data_input)
      
      # Reassign the original row names to the complete data
      rownames(new_data) <- rownames(data_input)[rownames(data_input) %in% rownames(new_data)]
      
      addMiss <- function(x, MV.rate, MNAR.ratio){
        
        rnorm <- rbinom <- na.omit <- NULL
        
        ini.seed = 233
        set.seed(ini.seed)
        seeds <- sample(1:10000,20)
        
        ## Add MNAR missing values according to protein abundance
        
        t_q <- stats::quantile(as.matrix(x), probs = MV.rate)  ## a th quantile of data matrix x
        
        set.seed(seeds[20])
        t_mat <- matrix(rnorm(n = nrow(x)*ncol(x), m = t_q, sd = 0.3),nrow = nrow(x),ncol = ncol(x))  ## threshold matrix
        
        set.seed(seeds[19])
        t_draw <- matrix(rbinom(nrow(x)*ncol(x), 1, MNAR.ratio), nrow = nrow(x),ncol = ncol(x))  ## Bernoulli draw matrix
        
        for (i in 1:nrow(x)) {
          for (j in 1:ncol(x)) {
            if (x[i,j] < t_mat[i,j] & t_draw[i,j] == 1) {
              x[i,j] <- NA
            }
          }
        }  ## add MNAR values
        
        ## sum(is.na(data1_1)) ## number of MNAR missing values
        
        if(sum(is.na(x)) > nrow(x)*ncol(x)*MV.rate) {return(x)}
        else{
          
          ## Add the rest missing values randomly
          
          x2 <- as.vector(as.matrix(x))
          ind1 <- which(is.na(x))
          ind2 <- seq(1:length(x2))
          ind2 <- ind2[!(ind2 %in% ind1)]
          set.seed(seeds[18])
          #ind <- which(x2 %in% sample(x2[!is.na(x2)], as.integer((nrow(x)*ncol(x)*MV.rate-sum(is.na(x))))))
          ind2 <- sample(ind2,as.integer((nrow(x)*ncol(x)*MV.rate-sum(is.na(x)))))
          ind <- sort(c(ind1,ind2))
          x2[ind] <- NA  ##  Add a*b MCAR values randomly
          ## sum(is.na(x2))  ## number of total missing values
          x2 <- matrix(x2, nrow = nrow(x),ncol = ncol(x))
          x2 <- as.data.frame(x2)
          colnames(x2) <- colnames(x)
          rownames(x2) <- rownames(x)
          return(x2)
        }
      }
      
      set.seed(123)
      options(digits = 5)
      
      # # p% MV
      # mv_p <- function(data, p){
      #   set.seed(300)
      #   addMiss(data, MV.rate = p, MNAR.ratio = 0)
      # }
      
      #10% MV - 0% MNAR
      mv10_0_new_data1 <- addMiss(new_data, 0.1,0)
      mv10_0_new_data<-cbind(row.names(new_data), mv10_0_new_data1)
      
      #10% MV - 20% MNAR
      mv10_20_new_data1 <- addMiss(new_data, 0.1,0.2)
      mv10_20_new_data<-cbind(row.names(new_data), mv10_20_new_data1)
      
      #10% MV - 50% MNAR
      mv10_50_new_data1 <- addMiss(new_data, 0.1,0.5)
      mv10_50_new_data<-cbind(row.names(new_data), mv10_50_new_data1)
      
      #10% MV - 80% MNAR
      mv10_80_new_data1 <- addMiss(new_data, 0.1,0.8)
      mv10_80_new_data<-cbind(row.names(new_data), mv10_80_new_data1)
      
      #20% MV - 0% MNAR
      mv20_0_new_data1 <- addMiss(new_data, 0.2,0)
      mv20_0_new_data<-cbind(row.names(new_data), mv20_0_new_data1)
      
      #20% MV - 20% MNAR
      mv20_20_new_data1 <- addMiss(new_data, 0.2,0.2)
      mv20_20_new_data<-cbind(row.names(new_data), mv20_20_new_data1)
      
      #20% MV - 50% MNAR
      mv20_50_new_data1 <- addMiss(new_data, 0.2,0.5)
      mv20_50_new_data<-cbind(row.names(new_data), mv20_50_new_data1)
      
      #20% MV - 80% MNAR
      mv20_80_new_data1 <- addMiss(new_data, 0.2,0.8)
      mv20_80_new_data<-cbind(row.names(new_data), mv20_80_new_data1)
      
      #30% MV - 0% MNAR
      mv30_0_new_data1 <- addMiss(new_data, 0.3,0)
      mv30_0_new_data<-cbind(row.names(new_data), mv30_0_new_data1)
      
      #30% MV - 20% MNAR
      mv30_20_new_data1 <- addMiss(new_data, 0.3,0.2)
      mv30_20_new_data<-cbind(row.names(new_data), mv30_20_new_data1)
      
      #30% MV - 50% MNAR
      mv30_50_new_data1 <- addMiss(new_data, 0.3,0.5)
      mv30_50_new_data<-cbind(row.names(new_data), mv30_50_new_data1)
      
      #30% MV - 80% MNAR
      mv30_80_new_data1 <- addMiss(new_data, 0.3,0.8)
      mv30_80_new_data<-cbind(row.names(new_data), mv30_80_new_data1)
      
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
      
      com_mv10_0_new_data1 <- complete_data_fn(mv10_0_new_data, groups)
      com_mv10_20_new_data1 <- complete_data_fn(mv10_20_new_data, groups)
      com_mv10_50_new_data1 <- complete_data_fn(mv10_50_new_data, groups)
      com_mv10_80_new_data1 <- complete_data_fn(mv10_80_new_data, groups)
      
      com_mv20_0_new_data1 <- complete_data_fn(mv20_0_new_data, groups)
      com_mv20_20_new_data1 <- complete_data_fn(mv20_20_new_data, groups)
      com_mv20_50_new_data1 <- complete_data_fn(mv20_50_new_data, groups)
      com_mv20_80_new_data1 <- complete_data_fn(mv20_80_new_data, groups)
      
      com_mv30_0_new_data1 <- complete_data_fn(mv30_0_new_data, groups)
      com_mv30_20_new_data1 <- complete_data_fn(mv30_20_new_data, groups)
      com_mv30_50_new_data1 <- complete_data_fn(mv30_50_new_data, groups)
      com_mv30_80_new_data1 <- complete_data_fn(mv30_80_new_data, groups)
      
      #Giving original name to the first column
      colnames(com_mv10_0_new_data1)[1] <- colnames(com_mv10_20_new_data1)[1]<-colnames(com_mv10_50_new_data1)[1]<-colnames(com_mv10_80_new_data1)[1]<-("rowid")
      colnames(com_mv20_0_new_data1)[1] <- colnames(com_mv20_20_new_data1)[1]<-colnames(com_mv20_50_new_data1)[1]<-colnames(com_mv20_80_new_data1)[1]<-("rowid")
      colnames(com_mv30_0_new_data1)[1] <- colnames(com_mv30_20_new_data1)[1]<-colnames(com_mv30_50_new_data1)[1]<-colnames(com_mv30_80_new_data1)[1]<-("rowid")
      
      #Removing the ID column and selecting remaining data
      col_to_rownames <- function(data){
        result <- data%>%
          tibble::column_to_rownames(var=colnames(data)[1])
      }
      
      com_mv10_0_new_data <- col_to_rownames(com_mv10_0_new_data1)
      com_mv10_20_new_data <- col_to_rownames(com_mv10_20_new_data1)
      com_mv10_50_new_data <- col_to_rownames(com_mv10_50_new_data1)
      com_mv10_80_new_data <- col_to_rownames(com_mv10_80_new_data1)
      
      com_mv20_0_new_data <- col_to_rownames(com_mv20_0_new_data1)
      com_mv20_20_new_data <- col_to_rownames(com_mv20_20_new_data1)
      com_mv20_50_new_data <- col_to_rownames(com_mv20_50_new_data1)
      com_mv20_80_new_data <- col_to_rownames(com_mv20_80_new_data1)
      
      com_mv30_0_new_data <- col_to_rownames(com_mv30_0_new_data1)
      com_mv30_20_new_data <- col_to_rownames(com_mv30_20_new_data1)
      com_mv30_50_new_data <- col_to_rownames(com_mv30_50_new_data1)
      com_mv30_80_new_data <- col_to_rownames(com_mv30_80_new_data1)
      
      
      #Grouping of dataframe
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
      group_data_mv10_0 <- grouping_data(com_mv10_0_new_data, groups)
      group_data_mv10_20 <- grouping_data(com_mv10_20_new_data, groups)
      group_data_mv10_50 <- grouping_data(com_mv10_50_new_data, groups)
      group_data_mv10_80 <- grouping_data(com_mv10_80_new_data, groups)
      
      group_data_mv20_0 <- grouping_data(com_mv20_0_new_data, groups)
      group_data_mv20_20 <- grouping_data(com_mv20_20_new_data, groups)
      group_data_mv20_50 <- grouping_data(com_mv20_50_new_data, groups)
      group_data_mv20_80 <- grouping_data(com_mv20_80_new_data, groups)
      
      group_data_mv30_0 <- grouping_data(com_mv30_0_new_data, groups)
      group_data_mv30_20 <- grouping_data(com_mv30_20_new_data, groups)
      group_data_mv30_50 <- grouping_data(com_mv30_50_new_data, groups)
      group_data_mv30_80 <- grouping_data(com_mv30_80_new_data, groups)
      
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
        newdata<- as.matrix(log2(dat))
        cycLoessNormed <- limma::normalizeCyclicLoess(newdata)
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
                               sampleLog2Median, init = "ls", psi = psi.huber,
                             na.action = stats::na.exclude, maxit = 20)
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
      #Imputation of normalized datasets
      #KNN imputation
      KNN_Imputation <- function(dat) {
        set.seed(123)
        
        # Preserve row names
        row_names <- rownames(dat)
        
        # Perform KNN imputation
        resultkNN <- VIM::kNN(dat, numFun = laeken::weightedMean, weightDist = TRUE,
                              imp_var = FALSE, k = 10)
        
        # Convert result to matrix to preserve row names
        result_matrix <- as.matrix(resultkNN)
        
        # Assign original row names to the result
        rownames(result_matrix) <- row_names
        
        return(result_matrix)
      }
      
      
      LLS_Imputation <- function(dat) {
        # Preserve row names
        row_names <- rownames(dat)
        
        # Perform LLS imputation
        resultLLS <- pcaMethods::llsImpute(as.matrix(dat), k = 2, correlation = "pearson", allVariables = TRUE)
        dataSet.imputed <- resultLLS@completeObs
        
        # Convert result to matrix if needed
        dataSet.imputed <- as.matrix(dataSet.imputed)
        
        # Assign original row names to the result
        rownames(dataSet.imputed) <- row_names
        
        return(dataSet.imputed)
      }
      
      SVD_Imputation <- function(dat) {
        # Preserve row names
        row_names <- rownames(dat)
        
        # Perform SVD imputation
        resultSVD <- pcaMethods::pca(dat, method = "svdImpute", nPcs = 2)
        dataSet.imputed <- resultSVD@completeObs
        
        # Convert result to matrix if needed
        dataSet.imputed <- as.matrix(dataSet.imputed)
        
        # Assign original row names to the result
        rownames(dataSet.imputed) <- row_names
        
        return(dataSet.imputed)
      }
      
      #vsn 
      vsn_mv <- function (data){
        result1 <- do.call("cbind", lapply(data, VSN_Norm))
        colnames(result1) <- colnames(new_data)
        result <- grouping_data(result1, groups)
        return(result)
      }
      vsn_group_data_mv10_0 <- vsn_mv(group_data_mv10_0)
      vsn_group_data_mv10_20 <- vsn_mv(group_data_mv10_20)
      vsn_group_data_mv10_50 <- vsn_mv(group_data_mv10_50)
      vsn_group_data_mv10_80 <- vsn_mv(group_data_mv10_80)
      
      vsn_group_data_mv20_0 <- vsn_mv(group_data_mv20_0)
      vsn_group_data_mv20_20 <- vsn_mv(group_data_mv20_20)
      vsn_group_data_mv20_50 <- vsn_mv(group_data_mv20_50)
      vsn_group_data_mv20_80 <- vsn_mv(group_data_mv20_80)
      
      vsn_group_data_mv30_0 <- vsn_mv(group_data_mv30_0)
      vsn_group_data_mv30_20 <- vsn_mv(group_data_mv30_20)
      vsn_group_data_mv30_50 <- vsn_mv(group_data_mv30_50)
      vsn_group_data_mv30_80 <- vsn_mv(group_data_mv30_80)
      
      #loess
      loess_mv <- function (data){
        result1 <- do.call("cbind", lapply(data, LOESS_Norm))
        colnames(result1) <- colnames(new_data)
        result <- grouping_data(result1, groups)
        return(result)
      }
      
      loess_group_data_mv10_0 <- loess_mv(group_data_mv10_0)
      loess_group_data_mv10_20 <- loess_mv(group_data_mv10_20)
      loess_group_data_mv10_50 <- loess_mv(group_data_mv10_50)
      loess_group_data_mv10_80 <- loess_mv(group_data_mv10_80)
      
      loess_group_data_mv20_0 <- loess_mv(group_data_mv20_0)
      loess_group_data_mv20_20 <- loess_mv(group_data_mv20_20)
      loess_group_data_mv20_50 <- loess_mv(group_data_mv20_50)
      loess_group_data_mv20_80 <- loess_mv(group_data_mv20_80)
      
      loess_group_data_mv30_0 <- loess_mv(group_data_mv30_0)
      loess_group_data_mv30_20 <- loess_mv(group_data_mv30_20)
      loess_group_data_mv30_50 <- loess_mv(group_data_mv30_50)
      loess_group_data_mv30_80 <- loess_mv(group_data_mv30_80)
      
      #rlr 
      rlr_mv <- function (data){
        result1 <- do.call("cbind", lapply(data, RLR_Norm))
        colnames(result1) <- colnames(new_data)
        result <- grouping_data(result1, groups)
        return(result)
      }
      
      rlr_group_data_mv10_0 <- rlr_mv(group_data_mv10_0)
      rlr_group_data_mv10_20 <- rlr_mv(group_data_mv10_20)
      rlr_group_data_mv10_50 <- rlr_mv(group_data_mv10_50)
      rlr_group_data_mv10_80 <- rlr_mv(group_data_mv10_80)
      
      rlr_group_data_mv20_0 <- rlr_mv(group_data_mv20_0)
      rlr_group_data_mv20_20 <- rlr_mv(group_data_mv20_20)
      rlr_group_data_mv20_50 <- rlr_mv(group_data_mv20_50)
      rlr_group_data_mv20_80 <- rlr_mv(group_data_mv20_80)
      
      rlr_group_data_mv30_0 <- rlr_mv(group_data_mv30_0)
      rlr_group_data_mv30_20 <- rlr_mv(group_data_mv30_20)
      rlr_group_data_mv30_50 <- rlr_mv(group_data_mv30_50)
      rlr_group_data_mv30_80 <- rlr_mv(group_data_mv30_80)
      
      #VSN and KNN 
      mv10_0_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv10_0, KNN_Imputation))
      mv10_20_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv10_20, KNN_Imputation))
      mv10_50_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv10_50, KNN_Imputation))
      mv10_80_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv10_80, KNN_Imputation))
      
      mv20_0_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv20_0, KNN_Imputation))
      mv20_20_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv20_20, KNN_Imputation))
      mv20_50_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv20_50, KNN_Imputation))
      mv20_80_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv20_80, KNN_Imputation))
      
      mv30_0_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv30_0, KNN_Imputation))
      mv30_20_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv30_20, KNN_Imputation))
      mv30_50_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv30_50, KNN_Imputation))
      mv30_80_vsn_knn_data <- do.call("cbind", lapply(vsn_group_data_mv30_80, KNN_Imputation))
      
      #VSN and LLS 
      mv10_0_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv10_0, LLS_Imputation))
      mv10_20_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv10_20, LLS_Imputation))
      mv10_50_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv10_50, LLS_Imputation))
      mv10_80_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv10_80, LLS_Imputation))
      
      mv20_0_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv20_0, LLS_Imputation))
      mv20_20_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv20_20, LLS_Imputation))
      mv20_50_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv20_50, LLS_Imputation))
      mv20_80_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv20_80, LLS_Imputation))
      
      mv30_0_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv30_0, LLS_Imputation))
      mv30_20_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv30_20, LLS_Imputation))
      mv30_50_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv30_50, LLS_Imputation))
      mv30_80_vsn_lls_data <- do.call("cbind", lapply(vsn_group_data_mv30_80, LLS_Imputation))
      
      #VSN and SVD 
      mv10_0_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv10_0, SVD_Imputation))
      mv10_20_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv10_20, SVD_Imputation))
      mv10_50_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv10_50, SVD_Imputation))
      mv10_80_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv10_80, SVD_Imputation))
      
      mv20_0_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv20_0, SVD_Imputation))
      mv20_20_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv20_20, SVD_Imputation))
      mv20_50_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv20_50, SVD_Imputation))
      mv20_80_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv20_80, SVD_Imputation))
      
      mv30_0_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv30_0, SVD_Imputation))
      mv30_20_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv30_20, SVD_Imputation))
      mv30_50_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv30_50, SVD_Imputation))
      mv30_80_vsn_svd_data <- do.call("cbind", lapply(vsn_group_data_mv30_80, SVD_Imputation))
      
      #LOESS and KNN 
      mv10_0_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv10_0, KNN_Imputation))
      mv10_20_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv10_20, KNN_Imputation))
      mv10_50_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv10_50, KNN_Imputation))
      mv10_80_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv10_80, KNN_Imputation))
      
      mv20_0_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv20_0, KNN_Imputation))
      mv20_20_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv20_20, KNN_Imputation))
      mv20_50_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv20_50, KNN_Imputation))
      mv20_80_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv20_80, KNN_Imputation))
      
      mv30_0_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv30_0, KNN_Imputation))
      mv30_20_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv30_20, KNN_Imputation))
      mv30_50_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv30_50, KNN_Imputation))
      mv30_80_loess_knn_data <- do.call("cbind", lapply(loess_group_data_mv30_80, KNN_Imputation))
      
      #LOESS and LLS 
      mv10_0_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv10_0, LLS_Imputation))
      mv10_20_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv10_20, LLS_Imputation))
      mv10_50_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv10_50, LLS_Imputation))
      mv10_80_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv10_80, LLS_Imputation))
      
      mv20_0_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv20_0, LLS_Imputation))
      mv20_20_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv20_20, LLS_Imputation))
      mv20_50_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv20_50, LLS_Imputation))
      mv20_80_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv20_80, LLS_Imputation))
      
      mv30_0_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv30_0, LLS_Imputation))
      mv30_20_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv30_20, LLS_Imputation))
      mv30_50_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv30_50, LLS_Imputation))
      mv30_80_loess_lls_data <- do.call("cbind", lapply(loess_group_data_mv30_80, LLS_Imputation))
      
      #LOESS and SVD 
      mv10_0_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv10_0, SVD_Imputation))
      mv10_20_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv10_20, SVD_Imputation))
      mv10_50_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv10_50, SVD_Imputation))
      mv10_80_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv10_80, SVD_Imputation))
      
      mv20_0_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv20_0, SVD_Imputation))
      mv20_20_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv20_20, SVD_Imputation))
      mv20_50_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv20_50, SVD_Imputation))
      mv20_80_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv20_80, SVD_Imputation))
      
      mv30_0_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv30_0, SVD_Imputation))
      mv30_20_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv30_20, SVD_Imputation))
      mv30_50_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv30_50, SVD_Imputation))
      mv30_80_loess_svd_data <- do.call("cbind", lapply(loess_group_data_mv30_80, SVD_Imputation))
      
      #RLR and KNN 
      mv10_0_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv10_0, KNN_Imputation))
      mv10_20_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv10_20, KNN_Imputation))
      mv10_50_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv10_50, KNN_Imputation))
      mv10_80_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv10_80, KNN_Imputation))
      
      mv20_0_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv20_0, KNN_Imputation))
      mv20_20_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv20_20, KNN_Imputation))
      mv20_50_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv20_50, KNN_Imputation))
      mv20_80_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv20_80, KNN_Imputation))
      
      mv30_0_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv30_0, KNN_Imputation))
      mv30_20_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv30_20, KNN_Imputation))
      mv30_50_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv30_50, KNN_Imputation))
      mv30_80_rlr_knn_data <- do.call("cbind", lapply(rlr_group_data_mv30_80, KNN_Imputation))
      
      #RLR and LLS 
      mv10_0_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv10_0, LLS_Imputation))
      mv10_20_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv10_20, LLS_Imputation))
      mv10_50_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv10_50, LLS_Imputation))
      mv10_80_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv10_80, LLS_Imputation))
      
      mv20_0_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv20_0, LLS_Imputation))
      mv20_20_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv20_20, LLS_Imputation))
      mv20_50_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv20_50, LLS_Imputation))
      mv20_80_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv20_80, LLS_Imputation))
      
      mv30_0_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv30_0, LLS_Imputation))
      mv30_20_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv30_20, LLS_Imputation))
      mv30_50_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv30_50, LLS_Imputation))
      mv30_80_rlr_lls_data <- do.call("cbind", lapply(rlr_group_data_mv30_80, LLS_Imputation))
      
      #RLR and SVD 
      mv10_0_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv10_0, SVD_Imputation))
      mv10_20_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv10_20, SVD_Imputation))
      mv10_50_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv10_50, SVD_Imputation))
      mv10_80_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv10_80, SVD_Imputation))
      
      mv20_0_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv20_0, SVD_Imputation))
      mv20_20_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv20_20, SVD_Imputation))
      mv20_50_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv20_50, SVD_Imputation))
      mv20_80_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv20_80, SVD_Imputation))
      
      mv30_0_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv30_0, SVD_Imputation))
      mv30_20_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv30_20, SVD_Imputation))
      mv30_50_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv30_50, SVD_Imputation))
      mv30_80_rlr_svd_data <- do.call("cbind", lapply(rlr_group_data_mv30_80, SVD_Imputation))
      
      new_data <- as.data.frame(new_data) %>%
        tibble::rownames_to_column(var="rowid")
      
      # Convert row names to a column named 'rowid'
      row_colname <- function (data){
        result <- as.data.frame(data)%>%
          tibble::rownames_to_column(var = "rowid")
      }
      mv10_0_vsn_knn_data <- row_colname(mv10_0_vsn_knn_data)
      mv10_20_vsn_knn_data <- row_colname(mv10_20_vsn_knn_data)
      mv10_50_vsn_knn_data <- row_colname(mv10_50_vsn_knn_data)
      mv10_80_vsn_knn_data <- row_colname(mv10_80_vsn_knn_data)
      
      mv10_0_vsn_lls_data <- row_colname(mv10_0_vsn_lls_data)
      mv10_20_vsn_lls_data <- row_colname(mv10_20_vsn_lls_data)
      mv10_50_vsn_lls_data <- row_colname(mv10_50_vsn_lls_data)
      mv10_80_vsn_lls_data <- row_colname(mv10_80_vsn_lls_data)
      
      mv10_0_vsn_svd_data <- row_colname(mv10_0_vsn_svd_data)
      mv10_20_vsn_svd_data <- row_colname(mv10_20_vsn_svd_data)
      mv10_50_vsn_svd_data <- row_colname(mv10_50_vsn_svd_data)
      mv10_80_vsn_svd_data <- row_colname(mv10_80_vsn_svd_data)
      
      mv20_0_vsn_knn_data <- row_colname(mv20_0_vsn_knn_data)
      mv20_20_vsn_knn_data <- row_colname(mv20_20_vsn_knn_data)
      mv20_50_vsn_knn_data <- row_colname(mv20_50_vsn_knn_data)
      mv20_80_vsn_knn_data <- row_colname(mv20_80_vsn_knn_data)
      
      mv20_0_vsn_lls_data <- row_colname(mv20_0_vsn_lls_data)
      mv20_20_vsn_lls_data <- row_colname(mv20_20_vsn_lls_data)
      mv20_50_vsn_lls_data <- row_colname(mv20_50_vsn_lls_data)
      mv20_80_vsn_lls_data <- row_colname(mv20_80_vsn_lls_data)
      
      mv20_0_vsn_svd_data <- row_colname(mv20_0_vsn_svd_data)
      mv20_20_vsn_svd_data <- row_colname(mv20_20_vsn_svd_data)
      mv20_50_vsn_svd_data <- row_colname(mv20_50_vsn_svd_data)
      mv20_80_vsn_svd_data <- row_colname(mv20_80_vsn_svd_data)
      
      mv30_0_vsn_knn_data <- row_colname(mv30_0_vsn_knn_data)
      mv30_20_vsn_knn_data <- row_colname(mv30_20_vsn_knn_data)
      mv30_50_vsn_knn_data <- row_colname(mv30_50_vsn_knn_data)
      mv30_80_vsn_knn_data <- row_colname(mv30_80_vsn_knn_data)
      
      mv30_0_vsn_lls_data <- row_colname(mv30_0_vsn_lls_data)
      mv30_20_vsn_lls_data <- row_colname(mv30_20_vsn_lls_data)
      mv30_50_vsn_lls_data <- row_colname(mv30_50_vsn_lls_data)
      mv30_80_vsn_lls_data <- row_colname(mv30_80_vsn_lls_data)
      
      mv30_0_vsn_svd_data <- row_colname(mv30_0_vsn_svd_data)
      mv30_20_vsn_svd_data <- row_colname(mv30_20_vsn_svd_data)
      mv30_50_vsn_svd_data <- row_colname(mv30_50_vsn_svd_data)
      mv30_80_vsn_svd_data <- row_colname(mv30_80_vsn_svd_data)
      
      mv10_0_loess_knn_data <- row_colname(mv10_0_loess_knn_data)
      mv10_20_loess_knn_data <- row_colname(mv10_20_loess_knn_data)
      mv10_50_loess_knn_data <- row_colname(mv10_50_loess_knn_data)
      mv10_80_loess_knn_data <- row_colname(mv10_80_loess_knn_data)
      
      mv10_0_loess_lls_data <- row_colname(mv10_0_loess_lls_data)
      mv10_20_loess_lls_data <- row_colname(mv10_20_loess_lls_data)
      mv10_50_loess_lls_data <- row_colname(mv10_50_loess_lls_data)
      mv10_80_loess_lls_data <- row_colname(mv10_80_loess_lls_data)
      
      mv10_0_loess_svd_data <- row_colname(mv10_0_loess_svd_data)
      mv10_20_loess_svd_data <- row_colname(mv10_20_loess_svd_data)
      mv10_50_loess_svd_data <- row_colname(mv10_50_loess_svd_data)
      mv10_80_loess_svd_data <- row_colname(mv10_80_loess_svd_data)
      
      mv20_0_loess_knn_data <- row_colname(mv20_0_loess_knn_data)
      mv20_20_loess_knn_data <- row_colname(mv20_20_loess_knn_data)
      mv20_50_loess_knn_data <- row_colname(mv20_50_loess_knn_data)
      mv20_80_loess_knn_data <- row_colname(mv20_80_loess_knn_data)
      
      mv20_0_loess_lls_data <- row_colname(mv20_0_loess_lls_data)
      mv20_20_loess_lls_data <- row_colname(mv20_20_loess_lls_data)
      mv20_50_loess_lls_data <- row_colname(mv20_50_loess_lls_data)
      mv20_80_loess_lls_data <- row_colname(mv20_80_loess_lls_data)
      
      mv20_0_loess_svd_data <- row_colname(mv20_0_loess_svd_data)
      mv20_20_loess_svd_data <- row_colname(mv20_20_loess_svd_data)
      mv20_50_loess_svd_data <- row_colname(mv20_50_loess_svd_data)
      mv20_80_loess_svd_data <- row_colname(mv20_80_loess_svd_data)
      
      mv30_0_loess_knn_data <- row_colname(mv30_0_loess_knn_data)
      mv30_20_loess_knn_data <- row_colname(mv30_20_loess_knn_data)
      mv30_50_loess_knn_data <- row_colname(mv30_50_loess_knn_data)
      mv30_80_loess_knn_data <- row_colname(mv30_80_loess_knn_data)
      
      mv30_0_loess_lls_data <- row_colname(mv30_0_loess_lls_data)
      mv30_20_loess_lls_data <- row_colname(mv30_20_loess_lls_data)
      mv30_50_loess_lls_data <- row_colname(mv30_50_loess_lls_data)
      mv30_80_loess_lls_data <- row_colname(mv30_80_loess_lls_data)
      
      mv30_0_loess_svd_data <- row_colname(mv30_0_loess_svd_data)
      mv30_20_loess_svd_data <- row_colname(mv30_20_loess_svd_data)
      mv30_50_loess_svd_data <- row_colname(mv30_50_loess_svd_data)
      mv30_80_loess_svd_data <- row_colname(mv30_80_loess_svd_data)
      
      mv10_0_rlr_knn_data <- row_colname(mv10_0_rlr_knn_data)
      mv10_20_rlr_knn_data <- row_colname(mv10_20_rlr_knn_data)
      mv10_50_rlr_knn_data <- row_colname(mv10_50_rlr_knn_data)
      mv10_80_rlr_knn_data <- row_colname(mv10_80_rlr_knn_data)
      
      mv10_0_rlr_lls_data <- row_colname(mv10_0_rlr_lls_data)
      mv10_20_rlr_lls_data <- row_colname(mv10_20_rlr_lls_data)
      mv10_50_rlr_lls_data <- row_colname(mv10_50_rlr_lls_data)
      mv10_80_rlr_lls_data <- row_colname(mv10_80_rlr_lls_data)
      
      mv10_0_rlr_svd_data <- row_colname(mv10_0_rlr_svd_data)
      mv10_20_rlr_svd_data <- row_colname(mv10_20_rlr_svd_data)
      mv10_50_rlr_svd_data <- row_colname(mv10_50_rlr_svd_data)
      mv10_80_rlr_svd_data <- row_colname(mv10_80_rlr_svd_data)
      
      mv20_0_rlr_knn_data <- row_colname(mv20_0_rlr_knn_data)
      mv20_20_rlr_knn_data <- row_colname(mv20_20_rlr_knn_data)
      mv20_50_rlr_knn_data <- row_colname(mv20_50_rlr_knn_data)
      mv20_80_rlr_knn_data <- row_colname(mv20_80_rlr_knn_data)
      
      mv20_0_rlr_lls_data <- row_colname(mv20_0_rlr_lls_data)
      mv20_20_rlr_lls_data <- row_colname(mv20_20_rlr_lls_data)
      mv20_50_rlr_lls_data <- row_colname(mv20_50_rlr_lls_data)
      mv20_80_rlr_lls_data <- row_colname(mv20_80_rlr_lls_data)
      
      mv20_0_rlr_svd_data <- row_colname(mv20_0_rlr_svd_data)
      mv20_20_rlr_svd_data <- row_colname(mv20_20_rlr_svd_data)
      mv20_50_rlr_svd_data <- row_colname(mv20_50_rlr_svd_data)
      mv20_80_rlr_svd_data <- row_colname(mv20_80_rlr_svd_data)
      
      mv30_0_rlr_knn_data <- row_colname(mv30_0_rlr_knn_data)
      mv30_20_rlr_knn_data <- row_colname(mv30_20_rlr_knn_data)
      mv30_50_rlr_knn_data <- row_colname(mv30_50_rlr_knn_data)
      mv30_80_rlr_knn_data <- row_colname(mv30_80_rlr_knn_data)
      
      mv30_0_rlr_lls_data <- row_colname(mv30_0_rlr_lls_data)
      mv30_20_rlr_lls_data <- row_colname(mv30_20_rlr_lls_data)
      mv30_50_rlr_lls_data <- row_colname(mv30_50_rlr_lls_data)
      mv30_80_rlr_lls_data <- row_colname(mv30_80_rlr_lls_data)
      
      mv30_0_rlr_svd_data <- row_colname(mv30_0_rlr_svd_data)
      mv30_20_rlr_svd_data <- row_colname(mv30_20_rlr_svd_data)
      mv30_50_rlr_svd_data <- row_colname(mv30_50_rlr_svd_data)
      mv30_80_rlr_svd_data <- row_colname(mv30_80_rlr_svd_data)
      
      #nrmse function
      nrmse_fun <- function(ximp, xmis){
        common_rows <- Reduce(intersect, list(ximp[[1]], xmis[[1]], new_data[[1]]))
        
        # Step 2: Subset dataframes to only include rows with common identifiers
        ximp_matched <- ximp[ximp[[1]] %in% common_rows, ]
        xmis_matched <- xmis[xmis[[1]] %in% common_rows, ]
        new_data_matched <- new_data[new_data[[1]] %in% common_rows, ]
        
        # Ensure data is sorted in the same order based on the first column
        ximp_matched <- ximp_matched[order(ximp_matched[[1]]), ]
        xmis_matched <- xmis_matched[order(xmis_matched[[1]]), ]
        new_data_matched <- new_data_matched[order(new_data_matched[[1]]), ]
        
        # Remove the first column (identifiers) as it's not used in NRMSE calculation
        ximp_matched <- ximp_matched[, -1]
        xmis_matched <- xmis_matched[, -1]
        new_data_numeric <- new_data_matched[, -1]
        
        # Step 3: Calculate NRMSE for the matched rows
        nrmse_value <- missForest::nrmse(2^as.matrix(ximp_matched), 
                                         as.matrix(xmis_matched), 
                                         as.matrix(new_data_numeric))
        
        # Print the NRMSE result
        print(nrmse_value)
      }
      
      #Calculating NRMSE by `nrmse` function of missForest R package
      #For 10mv-0mnar datasets
      mv10_0_1 <- nrmse_fun(mv10_0_vsn_knn_data, com_mv10_0_new_data1)
      mv10_0_2 <- nrmse_fun(mv10_0_vsn_lls_data, com_mv10_0_new_data1)
      mv10_0_3 <- nrmse_fun(mv10_0_vsn_svd_data, com_mv10_0_new_data1)
      
      mv10_0_4 <- nrmse_fun(mv10_0_loess_knn_data, com_mv10_0_new_data1)
      mv10_0_5 <- nrmse_fun(mv10_0_loess_lls_data, com_mv10_0_new_data1)
      mv10_0_6 <- nrmse_fun(mv10_0_loess_svd_data, com_mv10_0_new_data1)
      
      mv10_0_7 <- nrmse_fun(mv10_0_rlr_knn_data, com_mv10_0_new_data1)
      mv10_0_8 <- nrmse_fun(mv10_0_rlr_lls_data, com_mv10_0_new_data1)
      mv10_0_9 <- nrmse_fun(mv10_0_rlr_svd_data, com_mv10_0_new_data1)
      
      
      mv10_0_nrmse <- c(mv10_0_1, mv10_0_2, mv10_0_3,
                        mv10_0_4, mv10_0_5, mv10_0_6,
                        mv10_0_7, mv10_0_8, mv10_0_9)
      
      #For 10mv-20mnar datasets
      mv10_20_1 <- nrmse_fun(mv10_20_vsn_knn_data, com_mv10_20_new_data1)
      mv10_20_2 <- nrmse_fun(mv10_20_vsn_lls_data, com_mv10_20_new_data1)
      mv10_20_3 <- nrmse_fun(mv10_20_vsn_svd_data, com_mv10_20_new_data1)
      
      mv10_20_4 <- nrmse_fun(mv10_20_loess_knn_data, com_mv10_20_new_data1)
      mv10_20_5 <- nrmse_fun(mv10_20_loess_lls_data, com_mv10_20_new_data1)
      mv10_20_6 <- nrmse_fun(mv10_20_loess_svd_data, com_mv10_20_new_data1)
      
      mv10_20_7 <- nrmse_fun(mv10_20_rlr_knn_data, com_mv10_20_new_data1)
      mv10_20_8 <- nrmse_fun(mv10_20_rlr_lls_data, com_mv10_20_new_data1)
      mv10_20_9 <- nrmse_fun(mv10_20_rlr_svd_data, com_mv10_20_new_data1)
      
      
      mv10_20_nrmse <- c(mv10_20_1, mv10_20_2, mv10_20_3,
                         mv10_20_4, mv10_20_5, mv10_20_6,
                         mv10_20_7, mv10_20_8, mv10_20_9)
      
      #For 10mv-50mnar datasets
      mv10_50_1 <- nrmse_fun(mv10_50_vsn_knn_data, com_mv10_50_new_data1)
      mv10_50_2 <- nrmse_fun(mv10_50_vsn_lls_data, com_mv10_50_new_data1)
      mv10_50_3 <- nrmse_fun(mv10_50_vsn_svd_data, com_mv10_50_new_data1)
      
      mv10_50_4 <- nrmse_fun(mv10_50_loess_knn_data, com_mv10_50_new_data1)
      mv10_50_5 <- nrmse_fun(mv10_50_loess_lls_data, com_mv10_50_new_data1)
      mv10_50_6 <- nrmse_fun(mv10_50_loess_svd_data, com_mv10_50_new_data1)
      
      mv10_50_7 <- nrmse_fun(mv10_50_rlr_knn_data, com_mv10_50_new_data1)
      mv10_50_8 <- nrmse_fun(mv10_50_rlr_lls_data, com_mv10_50_new_data1)
      mv10_50_9 <- nrmse_fun(mv10_50_rlr_svd_data, com_mv10_50_new_data1)
      
      
      mv10_50_nrmse <- c(mv10_50_1, mv10_50_2, mv10_50_3,
                         mv10_50_4, mv10_50_5, mv10_50_6,
                         mv10_50_7, mv10_50_8, mv10_50_9)
      
      #For 10mv-80mnar datasets
      mv10_80_1 <- nrmse_fun(mv10_80_vsn_knn_data, com_mv10_80_new_data1)
      mv10_80_2 <- nrmse_fun(mv10_80_vsn_lls_data, com_mv10_80_new_data1)
      mv10_80_3 <- nrmse_fun(mv10_80_vsn_svd_data, com_mv10_80_new_data1)
      
      mv10_80_4 <- nrmse_fun(mv10_80_loess_knn_data, com_mv10_80_new_data1)
      mv10_80_5 <- nrmse_fun(mv10_80_loess_lls_data, com_mv10_80_new_data1)
      mv10_80_6 <- nrmse_fun(mv10_80_loess_svd_data, com_mv10_80_new_data1)
      
      mv10_80_7 <- nrmse_fun(mv10_80_rlr_knn_data, com_mv10_80_new_data1)
      mv10_80_8 <- nrmse_fun(mv10_80_rlr_lls_data, com_mv10_80_new_data1)
      mv10_80_9 <- nrmse_fun(mv10_80_rlr_svd_data, com_mv10_80_new_data1)
      
      
      mv10_80_nrmse <- c(mv10_80_1, mv10_80_2, mv10_80_3,
                         mv10_80_4, mv10_80_5, mv10_80_6,
                         mv10_80_7, mv10_80_8, mv10_80_9)
      
      #Calculating NRMSE by `nrmse` function of missForest R package
      #For 20mv-0mnar datasets
      mv20_0_1 <- nrmse_fun(mv20_0_vsn_knn_data, com_mv20_0_new_data1)
      mv20_0_2 <- nrmse_fun(mv20_0_vsn_lls_data, com_mv20_0_new_data1)
      mv20_0_3 <- nrmse_fun(mv20_0_vsn_svd_data, com_mv20_0_new_data1)
      
      mv20_0_4 <- nrmse_fun(mv20_0_loess_knn_data, com_mv20_0_new_data1)
      mv20_0_5 <- nrmse_fun(mv20_0_loess_lls_data, com_mv20_0_new_data1)
      mv20_0_6 <- nrmse_fun(mv20_0_loess_svd_data, com_mv20_0_new_data1)
      
      mv20_0_7 <- nrmse_fun(mv20_0_rlr_knn_data, com_mv20_0_new_data1)
      mv20_0_8 <- nrmse_fun(mv20_0_rlr_lls_data, com_mv20_0_new_data1)
      mv20_0_9 <- nrmse_fun(mv20_0_rlr_svd_data, com_mv20_0_new_data1)
      
      
      mv20_0_nrmse <- c(mv20_0_1, mv20_0_2, mv20_0_3,
                        mv20_0_4, mv20_0_5, mv20_0_6,
                        mv20_0_7, mv20_0_8, mv20_0_9)
      
      #For 20mv-20mnar datasets
      mv20_20_1 <- nrmse_fun(mv20_20_vsn_knn_data, com_mv20_20_new_data1)
      mv20_20_2 <- nrmse_fun(mv20_20_vsn_lls_data, com_mv20_20_new_data1)
      mv20_20_3 <- nrmse_fun(mv20_20_vsn_svd_data, com_mv20_20_new_data1)
      
      mv20_20_4 <- nrmse_fun(mv20_20_loess_knn_data, com_mv20_20_new_data1)
      mv20_20_5 <- nrmse_fun(mv20_20_loess_lls_data, com_mv20_20_new_data1)
      mv20_20_6 <- nrmse_fun(mv20_20_loess_svd_data, com_mv20_20_new_data1)
      
      mv20_20_7 <- nrmse_fun(mv20_20_rlr_knn_data, com_mv20_20_new_data1)
      mv20_20_8 <- nrmse_fun(mv20_20_rlr_lls_data, com_mv20_20_new_data1)
      mv20_20_9 <- nrmse_fun(mv20_20_rlr_svd_data, com_mv20_20_new_data1)
      
      
      mv20_20_nrmse <- c(mv20_20_1, mv20_20_2, mv20_20_3,
                         mv20_20_4, mv20_20_5, mv20_20_6,
                         mv20_20_7, mv20_20_8, mv20_20_9)
      
      #For 20mv-50mnar datasets
      mv20_50_1 <- nrmse_fun(mv20_50_vsn_knn_data, com_mv20_50_new_data1)
      mv20_50_2 <- nrmse_fun(mv20_50_vsn_lls_data, com_mv20_50_new_data1)
      mv20_50_3 <- nrmse_fun(mv20_50_vsn_svd_data, com_mv20_50_new_data1)
      
      mv20_50_4 <- nrmse_fun(mv20_50_loess_knn_data, com_mv20_50_new_data1)
      mv20_50_5 <- nrmse_fun(mv20_50_loess_lls_data, com_mv20_50_new_data1)
      mv20_50_6 <- nrmse_fun(mv20_50_loess_svd_data, com_mv20_50_new_data1)
      
      mv20_50_7 <- nrmse_fun(mv20_50_rlr_knn_data, com_mv20_50_new_data1)
      mv20_50_8 <- nrmse_fun(mv20_50_rlr_lls_data, com_mv20_50_new_data1)
      mv20_50_9 <- nrmse_fun(mv20_50_rlr_svd_data, com_mv20_50_new_data1)
      
      
      mv20_50_nrmse <- c(mv20_50_1, mv20_50_2, mv20_50_3,
                         mv20_50_4, mv20_50_5, mv20_50_6,
                         mv20_50_7, mv20_50_8, mv20_50_9)
      
      #For 20mv-80mnar datasets
      mv20_80_1 <- nrmse_fun(mv20_80_vsn_knn_data, com_mv20_80_new_data1)
      mv20_80_2 <- nrmse_fun(mv20_80_vsn_lls_data, com_mv20_80_new_data1)
      mv20_80_3 <- nrmse_fun(mv20_80_vsn_svd_data, com_mv20_80_new_data1)
      
      mv20_80_4 <- nrmse_fun(mv20_80_loess_knn_data, com_mv20_80_new_data1)
      mv20_80_5 <- nrmse_fun(mv20_80_loess_lls_data, com_mv20_80_new_data1)
      mv20_80_6 <- nrmse_fun(mv20_80_loess_svd_data, com_mv20_80_new_data1)
      
      mv20_80_7 <- nrmse_fun(mv20_80_rlr_knn_data, com_mv20_80_new_data1)
      mv20_80_8 <- nrmse_fun(mv20_80_rlr_lls_data, com_mv20_80_new_data1)
      mv20_80_9 <- nrmse_fun(mv20_80_rlr_svd_data, com_mv20_80_new_data1)
      
      
      mv20_80_nrmse <- c(mv20_80_1, mv20_80_2, mv20_80_3,
                         mv20_80_4, mv20_80_5, mv20_80_6,
                         mv20_80_7, mv20_80_8, mv20_80_9)
      
      #Calculating NRMSE by `nrmse` function of missForest R package
      #For 30mv-0mnar datasets
      mv30_0_1 <- nrmse_fun(mv30_0_vsn_knn_data, com_mv30_0_new_data1)
      mv30_0_2 <- nrmse_fun(mv30_0_vsn_lls_data, com_mv30_0_new_data1)
      mv30_0_3 <- nrmse_fun(mv30_0_vsn_svd_data, com_mv30_0_new_data1)
      
      mv30_0_4 <- nrmse_fun(mv30_0_loess_knn_data, com_mv30_0_new_data1)
      mv30_0_5 <- nrmse_fun(mv30_0_loess_lls_data, com_mv30_0_new_data1)
      mv30_0_6 <- nrmse_fun(mv30_0_loess_svd_data, com_mv30_0_new_data1)
      
      mv30_0_7 <- nrmse_fun(mv30_0_rlr_knn_data, com_mv30_0_new_data1)
      mv30_0_8 <- nrmse_fun(mv30_0_rlr_lls_data, com_mv30_0_new_data1)
      mv30_0_9 <- nrmse_fun(mv30_0_rlr_svd_data, com_mv30_0_new_data1)
      
      
      mv30_0_nrmse <- c(mv30_0_1, mv30_0_2, mv30_0_3,
                        mv30_0_4, mv30_0_5, mv30_0_6,
                        mv30_0_7, mv30_0_8, mv30_0_9)
      
      #For 30mv-20mnar datasets
      mv30_20_1 <- nrmse_fun(mv30_20_vsn_knn_data, com_mv30_20_new_data1)
      mv30_20_2 <- nrmse_fun(mv30_20_vsn_lls_data, com_mv30_20_new_data1)
      mv30_20_3 <- nrmse_fun(mv30_20_vsn_svd_data, com_mv30_20_new_data1)
      
      mv30_20_4 <- nrmse_fun(mv30_20_loess_knn_data, com_mv30_20_new_data1)
      mv30_20_5 <- nrmse_fun(mv30_20_loess_lls_data, com_mv30_20_new_data1)
      mv30_20_6 <- nrmse_fun(mv30_20_loess_svd_data, com_mv30_20_new_data1)
      
      mv30_20_7 <- nrmse_fun(mv30_20_rlr_knn_data, com_mv30_20_new_data1)
      mv30_20_8 <- nrmse_fun(mv30_20_rlr_lls_data, com_mv30_20_new_data1)
      mv30_20_9 <- nrmse_fun(mv30_20_rlr_svd_data, com_mv30_20_new_data1)
      
      
      mv30_20_nrmse <- c(mv30_20_1, mv30_20_2, mv30_20_3,
                         mv30_20_4, mv30_20_5, mv30_20_6,
                         mv30_20_7, mv30_20_8, mv30_20_9)
      
      #For 30mv-50mnar datasets
      mv30_50_1 <- nrmse_fun(mv30_50_vsn_knn_data, com_mv30_50_new_data1)
      mv30_50_2 <- nrmse_fun(mv30_50_vsn_lls_data, com_mv30_50_new_data1)
      mv30_50_3 <- nrmse_fun(mv30_50_vsn_svd_data, com_mv30_50_new_data1)
      
      mv30_50_4 <- nrmse_fun(mv30_50_loess_knn_data, com_mv30_50_new_data1)
      mv30_50_5 <- nrmse_fun(mv30_50_loess_lls_data, com_mv30_50_new_data1)
      mv30_50_6 <- nrmse_fun(mv30_50_loess_svd_data, com_mv30_50_new_data1)
      
      mv30_50_7 <- nrmse_fun(mv30_50_rlr_knn_data, com_mv30_50_new_data1)
      mv30_50_8 <- nrmse_fun(mv30_50_rlr_lls_data, com_mv30_50_new_data1)
      mv30_50_9 <- nrmse_fun(mv30_50_rlr_svd_data, com_mv30_50_new_data1)
      
      
      mv30_50_nrmse <- c(mv30_50_1, mv30_50_2, mv30_50_3,
                         mv30_50_4, mv30_50_5, mv30_50_6,
                         mv30_50_7, mv30_50_8, mv30_50_9)
      
      #For 30mv-80mnar datasets
      mv30_80_1 <- nrmse_fun(mv30_80_vsn_knn_data, com_mv30_80_new_data1)
      mv30_80_2 <- nrmse_fun(mv30_80_vsn_lls_data, com_mv30_80_new_data1)
      mv30_80_3 <- nrmse_fun(mv30_80_vsn_svd_data, com_mv30_80_new_data1)
      
      mv30_80_4 <- nrmse_fun(mv30_80_loess_knn_data, com_mv30_80_new_data1)
      mv30_80_5 <- nrmse_fun(mv30_80_loess_lls_data, com_mv30_80_new_data1)
      mv30_80_6 <- nrmse_fun(mv30_80_loess_svd_data, com_mv30_80_new_data1)
      
      mv30_80_7 <- nrmse_fun(mv30_80_rlr_knn_data, com_mv30_80_new_data1)
      mv30_80_8 <- nrmse_fun(mv30_80_rlr_lls_data, com_mv30_80_new_data1)
      mv30_80_9 <- nrmse_fun(mv30_80_rlr_svd_data, com_mv30_80_new_data1)
      
      
      mv30_80_nrmse <- c(mv30_80_1, mv30_80_2, mv30_80_3,
                         mv30_80_4, mv30_80_5, mv30_80_6,
                         mv30_80_7, mv30_80_8, mv30_80_9)
      
      mv10_nrmse_result <- cbind(as.data.frame(mv10_0_nrmse), as.data.frame(mv10_20_nrmse),as.data.frame(mv10_50_nrmse),as.data.frame(mv10_80_nrmse))
      mv20_nrmse_result <- cbind(as.data.frame(mv20_0_nrmse), as.data.frame(mv20_20_nrmse),as.data.frame(mv20_50_nrmse),as.data.frame(mv20_80_nrmse))
      mv30_nrmse_result <- cbind(as.data.frame(mv30_0_nrmse), as.data.frame(mv30_20_nrmse),as.data.frame(mv30_50_nrmse),as.data.frame(mv30_80_nrmse))
      
      nrmse_result <- cbind(Combination = c("vsn_knn", "vsn_lls", "vsn_svd",
                                            "loess_knn", "loess_lls", "loess_svd",
                                            "rlr_knn", "rlr_lls", "rlr_svd"), mv10_nrmse_result, mv20_nrmse_result, mv30_nrmse_result)
      
      result = list ("NRMSE result" = nrmse_result)
      
      return(result)
    })
  })
}

nrmse_result1 <- nrmse_group_mv_norm_imp(X50vs0_5, lfproQC::yeast_groups)

df_long3 <- tidyr::pivot_longer(nrmse_result1$`NRMSE result`, cols = c("mv10_0_nrmse", "mv10_20_nrmse", "mv10_50_nrmse", "mv10_80_nrmse",
                                                       "mv20_0_nrmse", "mv20_20_nrmse", "mv20_50_nrmse", "mv20_80_nrmse",
                                                       "mv30_0_nrmse", "mv30_20_nrmse", "mv30_50_nrmse", "mv30_80_nrmse",), names_to = "variable", values_to = "value")
# Convert the line plot into a bar plot
# Convert the line plot into a bar plot
ggplot2::ggplot(df_long3, aes(x = Combination, y = value, fill = variable, group = variable)) +
  geom_col(position = "dodge", width = 0.7) +  # Dodged bars, for side-by-side comparison
  labs(title = "NRMSE comparison for groupwise normalized MVs generated datasets", x = "Combination", y = "Value") +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(12)) +  # Custom color palette
  theme_minimal() +
  theme(legend.position = "right",             # Move the legend to the right
        legend.title = element_text(size = 12),# Make legend title bigger
        legend.text = element_text(size = 10), # Make legend text bigger
        axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for clarity

