library(nnls)
library(tidyverse)
#library(tidymodels)
library(readxl)
#library(openxlsx)



## Setup automatic nnls deconvolution process from excel sheet


#### TODO ####
#filter(rsquared > 0.7) #cannot do this as some are still false positive and need to replace formula later with 0 for nnls


##############

# Read the skyline output and convert some variables and mutate new
df <- readxl::read_excel(fileInput) |>
    mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) |>
    mutate(Area = as.numeric(Area)) |>
    mutate(Area = replace_na(Area, 0)) |>  # Replace missing values with 0
    #mutate(RatioQuanToQual = as.numeric(RatioQuanToQual)) |>
    #mutate(RatioQualToQuan = as.numeric(RatioQualToQuan)) |>
    mutate(C_part = str_extract(Molecule, "C\\d+"),  # Extract "C" followed by numbers
           Cl_part = str_extract(Molecule, "Cl\\d+"), # Extract "Cl" followed by numbers
           C_number = as.numeric(str_extract(C_part, "\\d+")), # Extract numeric values for sorting
           Cl_number = as.numeric(str_extract(Cl_part, "\\d+")),

           # Combine them into simpliefied annotation for PCAs
           PCA = str_c(C_part, Cl_part, sep = "")
    )

if (correctWithRS == "Yes" & any(df$`Molecule List` == "RS")){
df <- df |>
    group_by(`Replicate Name`) |>
    mutate(Area = Area / first(Area[`Molecule List`== "RS" & `Isotope Label Type` == "Quan"])) |>
    ungroup()
}

# Calculate the average blank value, should be based on each homologue
if (blankSubtraction == "Yes"){ #input$blankSubtraction == "Yes"

#creating a df_blank with average of all blanks for each CxCly group
df_blank <- df |>
    filter(`Sample Type` == "Blank") |>
    #filter(`Isotope Label Type` == "Quan") |>
    group_by(Molecule, `Molecule List`, `Isotope Label Type`) |>
    summarize(AverageBlank = mean(Area, na.rm = TRUE)) |>
    ungroup() |>
    filter(!`Molecule List` %in% c("IS", "RS", "VS")) #dont include internal standards

# joining with the original df to get the AverageBlank column for all samples and homologues
Skyline_output <- df |>
    full_join(df_blank) |>
    mutate(AverageBlank = replace_na(AverageBlank, 0)) |>
    mutate(Area =  case_when(`Sample Type` == "Unknown" ~ Area - AverageBlank, .default = Area)) |> #mutate only when Sample Type == Unknown otherwise also subtract standards and others with blank
    mutate(Area = ifelse(Area <0, 0, Area))

} else {
    Skyline_output <- df
}


##### PREPARE FOR DECONVOLUTION #######
CPs_standards <- Skyline_output |>
    filter(`Sample Type` == "Standard", #the stds are not blank corrected
           !`Molecule List` %in% c("IS", "RS", "VS"), # dont include IS, RS, VS
           `Isotope Label Type` == "Quan", # use only Quant ions
           !!sym(standardAnnoColumn) != "NA") |>
    group_by(!!sym(standardAnnoColumn), Molecule, `Molecule List`, C_number, Cl_number, PCA) |> #input$standardAnnoColumn, "Note" is default. !!sym(standardAnnoColumn) unquotes the string variable and converts it to a symbol that dplyr can understand within the group_by() function
    nest() |>
    mutate(models = map(data, ~lm(Area ~ `Analyte Concentration`, data = .x))) |>
    mutate(coef = map(models, coef)) |>
    mutate(Response_factor = map_dbl(models, ~ coef(.x)["`Analyte Concentration`"]))|> #get the slope
    mutate(intercept = map(coef, pluck("(Intercept)"))) |>
    mutate(rsquared = map(models, summary)) |> #first creat a data frame list with the model
    mutate(rsquared = map(rsquared, pluck("r.squared"))) |> # then pluck only the r.squared value
    select(-coef) |>  # remove coef variable since it has already been plucked
    unnest(c(Response_factor, intercept, rsquared)) |>  #removing the list type for these variables
    mutate(Response_factor = if_else(Response_factor < 0, 0, Response_factor)) |> # replace negative RF with 0
    mutate(rsquared = ifelse(is.nan(rsquared), 0, rsquared)) |>
    mutate(Response_factor = if_else(rsquared < removeRsquared, 0, Response_factor)) |> #keep RF only if rsquared is above removeRsquared input
    ungroup() |>
    #rename(Chain_length = C_part) |>
    group_by(!!sym(standardAnnoColumn), C_number) |> #grouping by the selected Note, #input$standardAnnoColumn
    mutate(Sum_response_factor_chainlength = sum(Response_factor, na.rm = TRUE)) |>
    ungroup()


# This is for mixtures, single chain stds will be added later

if(standardTypes =="Mixtures"){
# For SCCPs
CPs_standards_S <- CPs_standards |>
    filter(str_detect(!!sym(standardAnnoColumn), "S-")) |>
    mutate(Response_factor = if_else(C_number < 14, Response_factor, 0)) #Need to restrict to C10-C13, if <C10 then vSCCPs?

# For MCCPs
CPs_standards_M <- CPs_standards |>
    filter(str_detect(!!sym(standardAnnoColumn), "M-")) |>
    mutate(Response_factor = if_else(C_number >= 14 & C_number <= 17, Response_factor, 0))

# For LCCPs
CPs_standards_L <- CPs_standards |>
    filter(str_detect(!!sym(standardAnnoColumn), "L-")) |>
    mutate(Response_factor = if_else(C_number >= 18, Response_factor, 0)) #Need to restrict to C18-C30, if >C30 then vLCCPs?

# Combine groups, only including those that have data
CPs_standards_list <- list(CPs_standards_S, CPs_standards_M, CPs_standards_L)

# Filter out empty data frames before binding
CPs_standards <- bind_rows(Filter(function(x) nrow(x) > 0, CPs_standards_list))
}


CPs_samples <- Skyline_output |>  #-> Skyline_output()
    filter(`Sample Type` == "Unknown",
           !`Molecule List` %in% c("IS", "RS", "VS"), # dont include IS, RS, VS
           `Isotope Label Type` == "Quan") |>
    #mutate(Group = case_when(
    #       Chain_length <10 ~ "vS",  # Group vS <10
    #       Chain_length >= 10 & Chain_length <= 13 ~ "S",  # Group S for 10-13
    #       Chain_length >= 14 & Chain_length <= 17 ~ "M",  # Group M for 14-17
    #      Chain_length >= 18 & Chain_length <= 30 ~ "L",  # Group L for 18-30
    #       Chain_length >30 ~ "vL",  # Group vL >30
    #     .default = "Unknown" # Default case
    #)) |>
    group_by(`Replicate Name`) |>  # Group by Replicate Name and Group
    mutate(Relative_distribution = Area / sum(Area, na.rm = TRUE)) |>
    ungroup() |>  # Ungroup before dropping the Group column
    #select(-Group) |>  # Explicitly remove Group column
    select(`Replicate Name`, Molecule, Area, Relative_distribution) |>
    mutate(across(Relative_distribution, ~replace(., is.nan(.), 0)))  # Replace NaN with zero


CPs_standards_input <- CPs_standards |>
    select(Molecule, !!sym(standardAnnoColumn), Response_factor) |> #-> !!sym(input$standardAnnoColumn)
    pivot_wider(names_from = !!sym(standardAnnoColumn), values_from = "Response_factor") #-> !!sym(input$standardAnnoColumn)


CPs_samples_input <- CPs_samples |>
    select(Molecule, `Replicate Name`, Relative_distribution) |>
    pivot_wider(names_from = "Replicate Name", values_from = "Relative_distribution")


# This step ensures that all values are corresponding to the same molecule for std and sample
problems_inputs <- tryCatch({
    CPs_samples_input |> right_join(CPs_standards_input, by = "Molecule")
}, warning = function(w) {
    message("A warning occurred: ", conditionMessage(w))
    NULL
}, error = function(e) {
    message("An error occurred: ", conditionMessage(e))
    NULL
})
# Check if result is NULL to handle the case where an error or warning occurred
if (is.null(problems_inputs)) {
    message("The operation did not complete successfully.")
} else {
    message("No problems! Standard and samples corresponds to each other.")
}

############################################################################### DECONVOLUTION #############################################################################

# Ensure combined_matrix is correctly defined as a matrix prior to the deconvolution

# First populate combined_matrix with CPs_standards_input
combined_matrix <- CPs_standards_input  |>
    select(-Molecule) |>
    as.matrix()


# Ensure combined_sample is correctly defined with nested data frames prior to the deconvolution
combined_sample <- CPs_samples  |>
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



########################################################## Calculate the concentration in ng/uL ###############################################################

#Calculate the response of the standards

#Remove the replicate name to generate vectors:
deconv_coef_df_matrix<- deconv_coef_df |>
    column_to_rownames(var = "Replicate Name")
    #select(-`Replicate Name`)

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
Final_results <- merged_df |>
    mutate(ConcentrationDetailed = Relative_distribution * Concentration)

# View the result
print(Final_results)


################################################### FINAL RESULTS ####################################################################


# Perform operations to reorganize data
Final_results <- Final_results |>
    select(-Area, -Relative_distribution, -Calculated_RF, -Measured_Signal, -Concentration) |> # Remove unwanted columns
    pivot_wider(
        names_from = Replicate.Name,  # Use values from Replicate.Name as new column names
        values_from = ConcentrationDetailed  # Use values from ConcentrationDetailed to fill the new columns
    )

#reorganized_data <- Concentration  |>
#       unnest(c(data)) |>
#      distinct(`Replicate.Name`, `Molecule`, .keep_all = TRUE)  |>
#     pivot_wider(names_from = `Molecule`, values_from = `Concentration`)
#reorganized_data <- t(reorganized_data) #transpose

#Make the first row (replicate names) the column names
#colnames(reorganized_data) <- reorganized_data[1, ]
#Samples_Concentration <- reorganized_data[-1, ]
# Convert the result back to a data frame
#Samples_Concentration <- as.data.frame(Samples_Concentration)
#Samples_Concentration2 <- Samples_Concentration |>
#       mutate(Molecule = CPs_samples_input$Molecule)|>
#      relocate(Molecule, .before = everything())

### END: Deconvolution script


