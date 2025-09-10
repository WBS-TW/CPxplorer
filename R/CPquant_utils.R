# Various utilities and helper functions for CPquant


#########################################################################################################
#------------------------------------------ CPquant functions ------------------------------------------#
#########################################################################################################

### Perform deconvolution on a single data frame ###
perform_deconvolution <- function(df, combined_standard, CPs_standards_sum_RF) {

    df_matrix <- as.matrix(df)

    print(paste("df_matrix dimensions:", dim(df_matrix)))
    print(paste("combined_standard dimensions:", dim(combined_standard)))

    if (nrow(combined_standard) != nrow(df_matrix)) {
        stop("Dimensions of combined_standard and df are incompatible.")
    }

    # Reshape df_matrix if it has only one column or extract the first column if it has multiple
    if (ncol(df_matrix) == 1) {
    df_vector <- as.vector(df_matrix)
    } else {
        df_vector <- as.vector(df_matrix["Relative_Area"])  # Extract the first column for nnls
    }

    # Check for NA/NaN/Inf values in df_vector and combined_standard
    if (any(is.na(df_vector)) || any(is.nan(df_vector)) || any(is.infinite(df_vector))) {
        stop("df_vector contains NA/NaN/Inf values.")
    }

    if (any(is.na(combined_standard)) || any(is.nan(combined_standard)) || any(is.infinite(combined_standard))) {
        stop("combined_standard contains NA/NaN/Inf values.")
    }

    # Perform nnls
    # combined_standard is response factor and df_vector is Relative_Area
    deconv <- nnls::nnls(combined_standard, df_vector)

    # Extract deconvolution results
    deconv_coef <- deconv$x


    #Normalize the coefficients so they sum to 1
    if (sum(deconv_coef) > 0) {
     deconv_coef <- deconv_coef / sum(deconv_coef)
    }else(deconv_coef <- 0)

    # Calculate deconvolved values (which is the response factor) using matrix multiplication
    deconv_resolved <- combined_standard %*% deconv_coef

    # Calculate the sum of RF*frac
    sum_deconv_RF <- as.matrix(CPs_standards_sum_RF) %*% deconv_coef





    # Calculate the goodness of fit by coefficient of determination (R2)
    # Calculate the total sum of squares (SST)
    sst <- sum((df_vector - mean(df_vector))^2)
    # Calculate the residual sum of squares (SSR)
    ssr <- sum((df_vector - deconv_resolved/sum(deconv_resolved))^2)

    # Calculate R-squared coefficient of determination
    deconv_rsquared <- 1 - (ssr / sst)

    # Calculate Mean Squared Error (MSE)
    mse <- mean((df_vector - deconv_resolved/sum(deconv_resolved))^2)

    # Calculate Root Mean Squared Error (RMSE)
    rmse <- sqrt(mse)

    #Chiq-square test: ensure that values are positive for chi-square test
    # if (any(deconv_resolved < 0) || any(df_vector < 0)) {
    #     warning("Non-positive values found, skipping chi-square test")
    #     chisq_result <- NULL
    # } else {
    #     #adding a very small constant to avoid 0 values
    #     observed_corr <- df_vector + 1E-12
    #     predicted_corr <- deconv_resolved + 1E-12
    #     chisq_result <- chisq.test(x= observed_corr, p = predicted_corr/sum(predicted_corr), rescale.p = TRUE)
    # }


    # Kolmogorov-Smirnov Test
    #ks_result <- ks.test(deconv_resolved, df_vector)



    #combine results for output
    combined_standard_names <- colnames(combined_standard)

    names(deconv_coef) <- combined_standard_names


    return(list(
        sum_deconv_RF = sum_deconv_RF,
        deconv_coef = deconv_coef,
        deconv_resolved = deconv_resolved,
        deconv_rsquared = deconv_rsquared
        #chisq_result = chisq_result
    ))
}

################################################################################



