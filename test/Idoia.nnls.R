---
        title: "R Notebook"
output: html_notebook
editor_options: 
        chunk_output_type: console
---
        
        
        library(nnls)
library(tidyverse)
library(tidymodels)
library(readxl)
library(openxlsx)



## Setup automatic nnls deconvolution process from excel sheet


#### TODO ####
#filter(rsquared > 0.7) #cannot do this as some are still false positive and need to replace formula later with 0 for nnls 

# NEED TO ADD using case_when or elseif..
# if (sum(!is.na(filtered_dataA$Area)) == 0 || sum(!is.na(filtered_dataA$`Analyte Concentration`)) == 0) {
#       cat("No valid cases for fitting the model for molecule:", molecule_name, "\n")
# OR add a filter to remove those rsquared that are very low (perhaps below 0.7?)  

##############

# reading data from excel
Skyline_output <- read_excel("F:/LINKOPING/CP analysis/ALEX/Results/Skyline_Results_Alex.xlsx") |> 
#Skyline_output <- read_excel("F:/LINKOPING/Manuscripts/Skyline/Data/Fish_Anders/Skyline/AndersSamples_ResultsFromSkyline_Concentration.xlsx") |> 
#Skyline_output <- read_excel("F:/LINKOPING/Manuscripts/Skyline/Skyline/OrbitrapDust.xlsx") |> 
        mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) |> 
        mutate(Area = as.numeric(Area)) |> 
        mutate(Area = replace_na(Area, 0)) |> # Replace missing values in the Response_factor column with 0
        mutate(Area = replace_na(Area, 0)) #|>
        #mutate(RatioQuanToQual = as.numeric(RatioQuanToQual)) |> #--< convert to numeric #-->
        #mutate(RatioQualToQuan = as.numeric(RatioQualToQuan)) #--< convert to numeric #-->


#This is currently filtered for C10-C13 only to compare with the Perkons script. Will remove later or add as an argument in function
CPs_standards <- Skyline_output |> 
        filter(`Sample Type` == "Standard",
               Molecule != "IS",
               Molecule != "RS",
               `Isotope Label Type` == "Quan",
               Note != "NA") |> 
        group_by(Note, Molecule) |>
        mutate(rel_int = Area/sum(Area)) |> #why ius it needed? Maybe can be removed
        nest() |> 
        mutate(models = map(data, ~lm(Area ~ `Analyte Concentration`, data = .x))) |> 
        mutate(coef = map(models, coef)) |> 
        mutate(Response_factor = map(coef, pluck("`Analyte Concentration`"))) |> 
        mutate(intercept = map(coef, pluck("(Intercept)"))) |> 
        mutate(rsquared = map(models, summary)) |> 
        mutate(rsquared = map(rsquared, pluck("r.squared"))) |>
        select(-coef) |>  # remove coef variable since it has already been plucked
        unnest(c(Response_factor, intercept, rsquared)) |>  #removing the list type for these variables
        mutate(Response_factor = if_else(Response_factor < 0, 0, Response_factor)) |> # replace negative RF with 0
        mutate(rsquared = ifelse(is.nan(rsquared), 0, rsquared)) |> 
        mutate(Chain_length = paste0("C", str_extract(Molecule, "(?<=C)[^H]+"))) |> 
        #filter(Chain_length == "C10" | Chain_length == "C11" | Chain_length == "C12" | Chain_length == "C13") |> #this will be remove later or added as arg in fn
        ungroup() |> 
        group_by(Note, Chain_length) |> 
        mutate(Sum_response_factor_chainlength = sum(Response_factor, na.rm = TRUE)) |> 
        ungroup()

# For SCCPs
CPs_standardsS <- CPs_standards |> 
        filter(str_detect(Note, "S-")) |> 
        mutate(Response_factor = if_else(Chain_length %in% c("C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C30"), 0, Response_factor))

# For MCCPs
CPs_standardsM <- CPs_standards |> 
        filter(str_detect(Note, "M-")) |> 
        mutate(Response_factor = if_else(Chain_length %in% c("C10", "C11", "C12", "C13", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C30"), 0, Response_factor))

# For LCCPs
CPs_standardsL <- CPs_standards |> 
        filter(str_detect(Note, "L-")) |> 
        mutate(Response_factor = if_else(Chain_length %in% c("C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17"), 0, Response_factor))

# Combine groups, only including those that have data
CPs_standards_list <- list(CPs_standardsS, CPs_standardsM, CPs_standardsL)

# Filter out empty data frames before binding
CPs_standards <- bind_rows(Filter(function(x) nrow(x) > 0, CPs_standards_list))


CPs_samples <- Skyline_output |> 
        filter(`Sample Type` == "Unknown",
               Molecule != "IS",
               Molecule != "RS",
               `Isotope Label Type` == "Quan") |> 
        mutate(Chain_length = str_extract(Molecule, "(?<=C)[^H]+") |> as.numeric()) |>  # Extract number after "C"
        #mutate(Group = case_when(
                #Chain_length >= 10 & Chain_length <= 13 ~ "S",  # Group S for 10-13
                #Chain_length >= 14 & Chain_length <= 17 ~ "M",  # Group M for 14-17
                #Chain_length >= 18 & Chain_length <= 30 ~ "L",  # Group L for 18-30
                #TRUE ~ "Unknown"  # Default case
        #)) |> 
        #group_by(`Replicate Name`, Group) |>  # Group by Replicate Name and Group
        group_by(`Replicate Name`) |>
        mutate(Relative_distribution = Area / sum(Area, na.rm = TRUE)) |> 
        ungroup() |>  # Ungroup before dropping the Group column
        #select(-Group) |>  # Explicitly remove Group column
        select(`Replicate Name`, Molecule, Area, Relative_distribution) |>
        mutate(across(Relative_distribution, ~replace(., is.nan(.), 0)))  # Replace NaN with zero




CPs_standards_input <- CPs_standards |> 
        select(Molecule, Note, Response_factor) |> 
        pivot_wider(names_from = "Note", values_from = "Response_factor")

CPs_samples_input <- CPs_samples |> 
        select(Molecule, `Replicate Name`, Relative_distribution) |> 
        pivot_wider(names_from = "Replicate Name", values_from = "Relative_distribution")

# This step ensures that all values are corresponding to the same molecule for std and sample        
combined <- CPs_samples_input |> 
        right_join(CPs_standards_input, by = "Molecule")


############################################################################### DECONVOLUTION #############################################################################

# Ensure combined_matrix is correctly defined as a matrix prior to the deconvolution
combined_matrix <- CPs_standards_input |> 
        select(-Molecule) |> 
        as.matrix()

# Ensure combined_sample is correctly defined with nested data frames prior to the deconvolution
combined_sample <- CPs_samples |> 
        group_by(`Replicate Name`) |> 
        select(-Molecule, -Area) |> 
        nest() |> 
        ungroup()

# Function to perform deconvolution on a single data frame
perform_deconvolution <- function(df, combined_matrix) {
        df_matrix <- as.matrix(df)
        
        print(paste("df_matrix dimensions:", dim(df_matrix)))
        print(paste("combined_matrix dimensions:", dim(combined_matrix)))
        
        if (nrow(combined_matrix) != nrow(df_matrix)) {
                stop("Dimensions of combined_matrix and df are incompatible.")
        }
        
        # Reshape df_matrix if it has only one column or extract the first column if it has multiple
        if (ncol(df_matrix) == 1) { 
                df_vector <- as.vector(df_matrix)
        } else {
                df_vector <- as.vector(df_matrix[, 1])  # Extract the first column for nnls
        }
        
        # Check for NA/NaN/Inf values in df_vector and combined_matrix
        if (any(is.na(df_vector)) || any(is.nan(df_vector)) || any(is.infinite(df_vector))) {
                stop("df_vector contains NA/NaN/Inf values.")
        }
        
        if (any(is.na(combined_matrix)) || any(is.nan(combined_matrix)) || any(is.infinite(combined_matrix))) {
                stop("combined_matrix contains NA/NaN/Inf values.")
        }
        
        # Perform nnls
        deconv <- nnls(combined_matrix, df_vector)
        
        # Extract deconvolution results
        deconv_coef <- deconv$x
        
        # Normalize the coefficients so they sum to 100%
        if (sum(deconv_coef) > 0) {
                deconv_coef <- deconv_coef / sum(deconv_coef) * 100
        }
        
        # Calculate deconvolved resolved values
        deconv_resolved <- combined_matrix %*% deconv_coef
        
        # Ensure that values are positive for chi-square test
        if (any(deconv_resolved <= 0) || any(df_vector <= 0)) {
                warning("Non-positive values found, skipping chi-square test")
                chisq_result <- NULL
        } else {
                chisq_result <- chisq.test(deconv_resolved, p = df_vector, rescale.p = TRUE)
        }
        
        return(list(
                deconv_coef = deconv_coef,
                deconv_resolved = deconv_resolved,
                chisq_result = chisq_result
        ))
}


# Apply the perform_deconvolution function to each nested data frame
Deconvolution <- combined_sample |> 
        mutate(result = map(data, ~ perform_deconvolution(.x, combined_matrix)))

# Extract deconv_coef from results and create a new data frame
deconv_coef_df <- Deconvolution |> 
        mutate(deconv_coef = map(result, "deconv_coef")) |> 
        select(`Replicate Name`, deconv_coef) |> 
        unnest_wider(deconv_coef, names_sep = "_")

# View the result
print(deconv_coef_df)



########################################################## Calculate the concentration in ng/uL ###############################################################
#Calculate the response of the standards

#Remove the replicate name to generate vectors:
deconv_coef_df_matrix<- deconv_coef_df |> 
        select(-`Replicate Name`)

        # Initialize an empty vector to store the summed results
        sum_results <- numeric(nrow(deconv_coef_df_matrix))
        
        # Iterate through each row of deconv_coef_df_matrix
        for (i in 1:nrow(deconv_coef_df_matrix)) {
                
                # Extract row vector from deconv_coef_df_matrix
                deconv_coef_vector <- as.numeric(deconv_coef_df_matrix[i, ])
                
                # Make sure the combined_matrix is correctly prepared each time
                combined_matrix <- CPs_standards_input |> 
                        select(-Molecule) |> 
                        mutate(across(everything(), as.numeric)) |> 
                        as.matrix()
                
                # Perform element-wise multiplication
                result <- sweep(combined_matrix, 2, deconv_coef_vector, `*`)
                
                # Sum all the values in the matrix for this iteration
                sum_results[i] <- sum(result, na.rm = TRUE) / 100
        }
        
        # Create a data frame with summed values and replicate names
        sum_results_df <- data.frame(
                `Replicate Name` = deconv_coef_df$`Replicate Name`,  # Assuming "Replicate Name" column exists
                Calculated_RF = sum_results
        )
        
        # View the resulting data frame
        print(sum_results_df)
 
#Calculate the sum of the area in the samples 
        #Sum the area of the samples
               Area <- CPs_samples |> 
                group_by(`Replicate Name`) |> 
                summarize(Measured_Signal = sum(Area, na.rm = TRUE)) |> 
                rename(Replicate.Name = `Replicate Name`)
        
        
#Calculate the concentration in ng/uL
       Final_results <- full_join(sum_results_df, Area, by = "Replicate.Name") |> 
               mutate(Concentration=Measured_Signal/Calculated_RF)

       CPs_samples<- CPs_samples |> 
               rename(Replicate.Name = `Replicate Name`)
       # Ensure both data frames are correctly named
       # Merge CPs_samples with Final_results based on the common column 'Replicate.Name'
       merged_df <- CPs_samples |> 
               left_join(Final_results, by = "Replicate.Name")
       
       # Multiply the column Relative_distribution by the corresponding value in Final_results
       # Assuming the column in Final_results to multiply is named 'value_column' (replace with actual column name)
       merged_df <- merged_df |> 
               mutate(ConcentrationDetailed = Relative_distribution * Concentration)
       
       # View the result
       print(merged_df)
       
              
######################################################### SAVE RESULTS ###################################################################
       
# Specify the file path where you want to save the Excel file
 excel_file <- "F:/LINKOPING/Manuscripts/Skyline/Skyline/Samples_Concentration_Not100.xlsx"
       
# Write 'Samples_Concentration' to Excel
write.xlsx(merged_df, excel_file, rowNames = FALSE)
       
# Confirm message that the file has been saved
cat("Excel file saved:", excel_file, "\n")

################################################### FINAL RESULTS ####################################################################
#CPs_samples<-CPs_samples |> 
 #       rename(`Replicate.Name` = `Replicate Name`)

# Merge total_sums_df into CPs_samples based on Replicate Name
#Concentration <- CPs_samples  |> 
 #       left_join(total_sums_df, by = "Replicate.Name")  |> 
  #      mutate(Concentration = `Relative_distribution` * `Total.Sum`)

#print(Concentration)
#Concentration<-Concentration |> 
 #       group_by(Replicate.Name) |> 
  #      distinct( `Molecule`, Concentration) |> 
   #     nest()

# Perform operations to reorganize data
#reorganized_data <- Concentration  |> 
 #       unnest() |>  
  #      distinct(`Replicate.Name`, `Molecule`, .keep_all = TRUE)  |> 
   #     pivot_wider(names_from = `Molecule`, values_from = `Concentration`)
#reorganized_data <- t(reorganized_data) #transpose

#Make the first row (replicate names) the column names
#colnames(reorganized_data) <- reorganized_data[1, ]
#Samples_Concentration <- reorganized_data[-1, ]
# Convert the result back to a data frame
#Samples_Concentration <- as.data.frame(Samples_Concentration)
#Samples_Concentration<- Samples_Concentration |> 
 #       mutate(Molecule = CPs_samples_input$Molecule)|> 
  #      relocate(Molecule, .before = everything()) 



















