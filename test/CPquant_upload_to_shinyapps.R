
library(readxl)
library(nnls)
library(dplyr)
library(tibble)
library(tidyr)
library(shiny)
library(DT)
library(plotly)
library(purrr)



options(shiny.maxRequestSize = 500 * 1024^2)
# UI
# UI
ui <- shiny::navbarPage("Quantification by deconvolution from Skyline output",
                        shiny::tabPanel("Quantification Inputs",
                                        shiny::fluidPage(shiny::sidebarLayout(
                                            shiny::sidebarPanel(
                                                shiny::fileInput("fileInput", "Import excel file from Skyline",
                                                                 accept = c('xlsx')),
                                                shiny::radioButtons("blankSubtraction",
                                                                    label = "Subtraction with blank?",
                                                                    choices = c("Yes", "No"), selected = "No"),
                                                shiny::radioButtons("correctWithRS", label = "Correct with RS?",
                                                                    choices = c("Yes", "No"), selected = "No"),
                                                shiny::radioButtons("standardTypes", label = "Types of standards",
                                                                    choices = c("Mixtures"), selected = "Mixtures"), #will add single chain std later
                                                shiny::actionButton('go', 'Proceed', width = "100%"),
                                                shiny::uiOutput("defineVariables")
                                            ),
                                            shiny::mainPanel(
                                                DT::DTOutput("table1")
                                            )
                                        )
                                        )),
                        shiny::tabPanel(
                            "Input summary",
                            shiny::mainPanel(
                                shiny::uiOutput("inputSummary")
                            )
                        ),

                        shiny::tabPanel(
                            "Quantification summary",
                            # shiny::sidebarLayout(
                            #         shiny::sidebarPanel(),
                            #         shiny::mainPanel(
                            #                 DT::DTOutput("quantTable")
                            #         )
                            # )
                            shiny::fluidPage(DT::DTOutput("quantTable"))
                        ),
                        shiny::tabPanel(
                            "QA/QC",
                            shiny::mainPanel(
                                DT::DTOutput("table2"),   # First table output (Skyline recovery data)
                                br(),                     # Optional line break to add space between the tables
                                DT::DTOutput("LOD")       # Second table output (LOD table)
                            )
                        ),
                        shiny::tabPanel(
                            "Instructions",
                            shiny::sidebarLayout(
                                shiny::sidebarPanel(shiny::h3("Manual"),
                                                    width = 3),
                                shiny::mainPanel(
                                    shiny::includeMarkdown("R/instructions_CPquant.md")
                                )
                            )
                        )

)


################################################################################
server <- function(input, output, session) {

    # Create a function to process each row

    Skyline_output <- reactive({
        req(input$fileInput) #requires that the input is available
        df <- readxl::read_excel(input$fileInput$datapath) |>
            dplyr::mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) |>
            dplyr::mutate(Area = as.numeric(Area)) |>
            dplyr::mutate(Area = replace_na(Area, 0)) |>  # Replace missing values with 0
            #dplyr::mutate(RatioQuanToQual = as.numeric(RatioQuanToQual)) |>
            #dplyr::mutate(RatioQualToQuan = as.numeric(RatioQualToQuan)) |>
            dplyr::mutate(C_part = stringr::str_extract(Molecule, "C\\d+"),  # Extract "C" followed by numbers
                          Cl_part = stringr::str_extract(Molecule, "Cl\\d+"), # Extract "Cl" followed by numbers
                          C_number = as.numeric(stringr::str_extract(C_part, "\\d+")), # Extract numeric values for sorting
                          Cl_number = as.numeric(stringr::str_extract(Cl_part, "\\d+")),
                          PCA = stringr::str_c(C_part, Cl_part, sep = "") # Combine them into simplified annotation for PCAs
            )



        #  Normalize data based on 'Correct with RS' input
        if (input$correctWithRS == "Yes" & any(df$`Molecule List` == "RS")){
            df <- df |>
                dplyr::group_by(`Replicate Name`) |>
                dplyr::mutate(Area = Area / first(Area[`Molecule List`== "RS" & `Isotope Label Type` == "Quan"])) |>
                dplyr::ungroup()
        }


        # Calculate the average blank value, should be based on each homologue
        if (input$blankSubtraction == "Yes"){

            #creating a df_blank with average of all blanks for each CxCly group
            df_blank <- df |>
                dplyr::filter(`Sample Type` == "Blank") |>
                #dplyr::filter(`Isotope Label Type` == "Quan") |>
                dplyr::group_by(Molecule, `Molecule List`, `Isotope Label Type`) |>
                dplyr::summarize(AverageBlank = mean(Area, na.rm = TRUE)) |>
                dplyr::ungroup() |>
                dplyr::filter(!`Molecule List` %in% c("IS", "RS", "VS")) #dont include internal standards

            # joining with the original df to get the AverageBlank column for all samples and homologues
            df <- df |>
                dplyr::full_join(df_blank) |>
                dplyr::mutate(AverageBlank = tidyr::replace_na(AverageBlank, 0)) |>
                dplyr::mutate(Area =  dplyr::case_when(`Sample Type` == "Unknown" ~ Area - AverageBlank, .default = Area)) |> #subtract only when Sample Type == Unknown otherwise also subtract standards and others with blank
                dplyr::mutate(Area = ifelse(Area <0, 0, Area)) #remove negative areas after blank subtraction

        } else {
            df
        }
    })



    ###START: Define UI components
    output$defineVariables <- shiny::renderUI({



        # Create the UI components
        shiny::fluidRow(
            shiny::h4("Define variables"),
            shiny::tags$br(),
            shiny::column(
                6,
                shiny::varSelectInput(
                    inputId = "standardAnnoColumn", #select which variable to use to define standards
                    label = "Variable for annotating standards",
                    data = Skyline_output(),
                    selected = "Batch Name"
                )
            ),
            # shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
            #shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
            #shiny::column(
            #       6,
            #      shiny::selectInput(
            #             inputId = "blanks", #select which variable to use to define standards
            #            label = "Define which samples are blanks",
            #           choices = unique(Skyline_output()$`Replicate Name`),
            #          multiple = TRUE
            # )
            #),
            shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
            shiny::column(
                6,
                shiny::selectInput(
                    inputId = "removeSamples", #select if some samples will be removed from quantification
                    label = 'Samples to remove from quantification?',
                    choices = unique(Skyline_output()$`Replicate Name`),
                    selected = NULL,
                    multiple = TRUE
                )
            ),
            # shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
            # shiny::column(
            # 6,
            #shiny::sliderInput(
            #       inputId = "removeAreas", #remove low peak areas
            #      label = "Keep absolute peak areas above this threshold (0 means keep everything)",
            #     min = min(Skyline_output()$Area),
            #    max = max(Skyline_output()$Area),
            #   value = 0,
            #  step = 10
            # )
            #),
            shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
            shiny::column(
                6,
                shiny::sliderInput(
                    inputId = "removeRsquared", #keep only Molecule above this rsquared, zero means keep everything
                    label = 'Keep the the calibration curves that show rsquared above this threshold (0 means keep everything)',
                    min = 0,
                    max = 1,
                    value = 0.80,
                    step = 0.05
                )
            )
        )
    })

    ### END: Define input variables

    # Set reactive values from user input
    # removeAreas <- eventReactive(input$go, {as.numeric(input$removeAreas)})
    removeRsquared <- shiny::eventReactive(input$go, {as.numeric(input$removeRsquared)})
    standardAnnoColumn <- shiny::eventReactive(input$go, {as.character(input$standardAnnoColumn)})
    removeSamples <- shiny::eventReactive(input$go, {as.character(input$removeSamples)})


    #Render raw table
    output$table1 <- DT::renderDT({
        DT::datatable(Skyline_output(),
                      options = list(
                          paging = TRUE,
                          pageLength = 50
                      )
        )
    })


    # Render summary statistics and plots of raw input BEFORE quantification
    output$inputSummary <- shiny::renderUI({
        shiny::fluidRow(
            shiny::column(6, plotly::plotlyOutput("plotSummary1", width = "30vw")),
            shiny::column(6, plotly::plotlyOutput("CalibrationSCCPs", width = "30vw")),
            shiny::column(6, plotly::plotlyOutput("CalibrationMCCPs", width = "30vw")),
            shiny::column(6, plotly::plotlyOutput("CalibrationLCCPs", width = "30vw"))
        )
    })


    output$plotSummary1 <- plotly::renderPlotly({
        data <- Skyline_output() |>
            dplyr::filter(`Isotope Label Type` == "Quan") |>
            dplyr::mutate(OrderedMolecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)]))) # Create a composite ordering factor

        plotly::plot_ly(
            data,
            x = ~ OrderedMolecule,
            y = ~ Area,
            color = ~ `Sample Type`,
            type = "box"
        )
    })


    # add load_plots() here



    ##### START: Deconvolution script

    shiny::observeEvent(input$go, {

        if(!is.null(removeSamples()) && length(removeSamples()) > 0){
            Skyline_output_filt <- Skyline_output() |>
                dplyr::filter(!`Replicate Name` %in% removeSamples())
        } else{
            Skyline_output_filt <- Skyline_output()
        }


        ##### PREPARE FOR DECONVOLUTION #######

        CPs_standards <- Skyline_output_filt |>
            dplyr::filter(`Sample Type` == "Standard", #make sure the stds are not blank corrected
                          !`Molecule List` %in% c("IS", "RS", "VS"), # dont include IS, RS, VS
                          `Isotope Label Type` == "Quan", # use only Quant ions
                          !!dplyr::sym(standardAnnoColumn()) != "NA") |>
            dplyr::group_by(!!dplyr::sym(standardAnnoColumn()), Molecule, `Molecule List`, C_number, Cl_number, PCA) |> #"Note" is default. !!sym(standardAnnoColumn) unquotes the string variable and converts it to a symbol that dplyr can understand within the group_by() function
            tidyr::nest() |>
            dplyr::mutate(models = purrr::map(data, ~lm(Area ~ `Analyte Concentration`, data = .x))) |>
            dplyr::mutate(coef = purrr::map(models, coef)) |>
            dplyr::mutate(Response_factor = purrr::map_dbl(models, ~ coef(.x)["`Analyte Concentration`"]))|> #get the slope
            dplyr::mutate(intercept = purrr::map(coef, purrr::pluck("(Intercept)"))) |>
            dplyr::mutate(rsquared = purrr::map(models, summary)) |> #first creat a data frame list with the model
            dplyr::mutate(rsquared = purrr::map(rsquared, purrr::pluck("r.squared"))) |> # then pluck only the r.squared value
            dplyr::select(-coef) |>  # remove coef variable since it has already been plucked
            tidyr::unnest(c(Response_factor, intercept, rsquared)) |>  #removing the list type for these variables
            dplyr::mutate(Response_factor = if_else(Response_factor < 0, 0, Response_factor)) |> # replace negative RF with 0
            dplyr::mutate(rsquared = ifelse(is.nan(rsquared), 0, rsquared)) |>
            dplyr::mutate(Response_factor = if_else(rsquared < removeRsquared(), 0, Response_factor)) |> #keep RF only if rsquared is above removeRsquared input
            dplyr::ungroup() |>
            dplyr::group_by(!!dplyr::sym(standardAnnoColumn()), C_number) |> #grouping by the selected Note, #input$standardAnnoColumn
            dplyr::mutate(Sum_response_factor_chainlength = sum(Response_factor, na.rm = TRUE)) |>
            dplyr::ungroup()



        # This is for mixtures, single chain stds will be added later

        if(input$standardTypes == "Mixtures"){
            # For SCCPs
            CPs_standards_S <- CPs_standards |>
                dplyr::filter(stringr::str_detect(!!dplyr::sym(standardAnnoColumn()), "S-")) |>
                dplyr::mutate(Response_factor = if_else(C_number < 14, Response_factor, 0)) #Need to restrict to C10-C13, if <C10 then vSCCPs?

            # For MCCPs
            CPs_standards_M <- CPs_standards |>
                dplyr::filter(stringr::str_detect(!!dplyr::sym(standardAnnoColumn()), "M-")) |>
                dplyr::mutate(Response_factor = if_else(C_number >= 14 & C_number <= 17, Response_factor, 0))

            # For LCCPs
            CPs_standards_L <- CPs_standards |>
                dplyr::filter(stringr::str_detect(!!dplyr::sym(standardAnnoColumn()), "L-")) |>
                dplyr::mutate(Response_factor = if_else(C_number >= 18, Response_factor, 0)) #Need to restrict to C18-C30, if >C30 then vLCCPs?

            # Combine groups, only including those that have data
            CPs_standards_list <- list(CPs_standards_S, CPs_standards_M, CPs_standards_L)

            # Filter out empty data frames before binding
            CPs_standards <- dplyr::bind_rows(Filter(function(x) nrow(x) > 0, CPs_standards_list))
        }



        CPs_samples <- Skyline_output_filt |>  #-> Skyline_output()
            dplyr::filter(`Sample Type` == "Unknown",
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
            dplyr::group_by(`Replicate Name`) |>  # Group by Replicate Name and Group
            dplyr::mutate(Relative_distribution = Area / sum(Area, na.rm = TRUE)) |>
            dplyr::ungroup() |>  # Ungroup before dropping the Group column
            #select(-Group) |>  # Explicitly remove Group column
            dplyr::select(`Replicate Name`, Molecule, Area, Relative_distribution) |>
            dplyr::mutate(dplyr::across(Relative_distribution, ~replace(., is.nan(.), 0)))  # Replace NaN with zero


        CPs_standards_input <- CPs_standards |>
            dplyr::select(Molecule, !!dplyr::sym(standardAnnoColumn()), Response_factor) |> #-> !!sym(input$standardAnnoColumn)
            tidyr::pivot_wider(names_from = !!dplyr::sym(standardAnnoColumn()), values_from = "Response_factor") #-> !!sym(input$standardAnnoColumn)


        CPs_samples_input <- CPs_samples |>
            dplyr::select(Molecule, `Replicate Name`, Relative_distribution) |>
            tidyr::pivot_wider(names_from = "Replicate Name", values_from = "Relative_distribution")


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
            dplyr::select(-Molecule) |>
            as.matrix()


        # Ensure combined_sample is correctly defined with nested data frames prior to the deconvolution
        combined_sample <- CPs_samples  |>
            dplyr::group_by(`Replicate Name`) |>
            dplyr::select(-Molecule, -Area) |>
            tidyr::nest() |>
            dplyr::ungroup()


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
            deconv <- nnls::nnls(combined_matrix, df_vector)

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
            dplyr::mutate(result = purrr::map(data, ~ perform_deconvolution(.x, combined_matrix)))

        # Extract deconv_coef from results and create a new data frame
        deconv_coef_df <- Deconvolution |>
            dplyr::mutate(deconv_coef = purrr::map(result, "deconv_coef")) |>
            dplyr::select(`Replicate Name`, deconv_coef) |>
            tidyr::unnest_wider(deconv_coef, names_sep = "_")


        ########################################################## Calculate the concentration in ng/uL ###############################################################

        #Calculate the response of the standards

        #Remove the replicate name to generate vectors:
        deconv_coef_df_matrix<- deconv_coef_df |>
            tibble::column_to_rownames(var = "Replicate Name")
        #select(-`Replicate Name`)

        # Initialize an empty vector to store the summed results
        sum_results <- numeric(nrow(deconv_coef_df_matrix))

        # Iterate through each row of deconv_coef_df_matrix
        for (i in 1:nrow(deconv_coef_df_matrix)) {

            # Extract row vector from deconv_coef_df_matrix
            deconv_coef_vector <- as.numeric(deconv_coef_df_matrix[i, ])

            # Make sure the combined_matrix is correctly prepared each time
            combined_matrix <- CPs_standards_input |>
                dplyr::select(-Molecule) |>
                dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
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
            dplyr::group_by(`Replicate Name`) |>
            dplyr::summarize(Measured_Signal = sum(Area, na.rm = TRUE)) |>
            dplyr::rename(Replicate.Name = `Replicate Name`)


        #Calculate the concentration in ng/uL
        Final_results <- dplyr::full_join(sum_results_df, Area, by = "Replicate.Name") |>
            dplyr::mutate(Concentration=Measured_Signal/Calculated_RF)

        CPs_samples<- CPs_samples |>
            dplyr::rename(Replicate.Name = `Replicate Name`)
        # Ensure both data frames are correctly named
        # Merge CPs_samples with Final_results based on the common column 'Replicate.Name'
        merged_df <- CPs_samples |>
            dplyr::left_join(Final_results, by = "Replicate.Name")

        # Multiply the column Relative_distribution by the corresponding value in Final_results
        # Assuming the column in Final_results to multiply is named 'value_column' (replace with actual column name)
        Final_results <- merged_df |>
            dplyr::mutate(ConcentrationDetailed = Relative_distribution * Concentration)

        # View the result
        print(Final_results)


        ################################################### FINAL RESULTS ####################################################################


        # Perform operations to reorganize data
        Final_results <- Final_results |>
            dplyr::select(-Area, -Relative_distribution, -Calculated_RF, -Measured_Signal, -Concentration) |> # Remove unwanted columns
            tidyr::pivot_wider(
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


        # Render table
        output$quantTable <- DT::renderDT({
            DT::datatable(Final_results,
                          filter = "top", extensions = c("Buttons", "Scroller"),
                          options = list(scrollY = 650,
                                         scrollX = 500,
                                         deferRender = TRUE,
                                         scroller = TRUE,
                                         buttons = list(list(extend = "excel", filename = "Samples_concentration", title = NULL,
                                                             exportOptions = list(
                                                                 modifier = list(page = "all")
                                                             )),
                                                        list(extend = "csv", filename = "Samples_concentration", title = NULL,
                                                             exportOptions = list(
                                                                 modifier = list(page = "all")
                                                             )),
                                                        list(extend = "colvis", targets = 0, visible = FALSE)),
                                         dom = "lBfrtip",
                                         fixedColumns = TRUE),
                          rownames = FALSE)

        })
    })


    ###########################################################RECOVERY#######################################################


    # Define the QA/QC reactive block
    Skyline_recovery <- reactive({
        # Ensure data is available
        df <- Skyline_output_filt  # Use the reactive data source
        req(df)  # Make sure df is not NULL

        # Prepare RECOVERY Data
        RECOVERY <- df |>
            dplyr::filter(`Isotope Label Type` == "Quan") |>
            tidyr::pivot_wider(id_cols = c(`Replicate Name`, `Sample Type`),
                               names_from = Molecule,
                               values_from = Area) |>  # Spread IS and RS into columns
            as.data.frame()  # Convert to data frame for consistency

        # Check if IS and RS columns exist before trying to mutate them
        if ("IS" %in% colnames(RECOVERY) & "RS" %in% colnames(RECOVERY)) {
            RECOVERY <- RECOVERY |>
                dplyr::mutate(across(c(IS, RS), ~replace_na(.x, 0)))  # Replace NAs with 0 for IS and RS
        } else {
            stop("Required columns 'IS' or 'RS' are missing from the data.")
        }

        # Calculate AverageRatio for Quality Control
        dfR <- RECOVERY |>
            dplyr::filter(`Sample Type` == "Quality Control") |>
            dplyr::mutate(RatioStd = IS / RS) |>
            dplyr::summarize(AverageRatio = mean(RatioStd, na.rm = TRUE))

        # Ensure dfR is valid
        if (nrow(dfR) == 0) {
            stop("No Quality Control samples found.")
        }

        # RECOVERY Calculation
        RECOVERY <- RECOVERY |>
            dplyr::filter(`Sample Type` %in% c("Unknown", "Blank")) |>
            dplyr::mutate(RatioSample = IS / RS) |>
            dplyr::mutate(Recovery = RatioSample / as.numeric(dfR$AverageRatio)) |>
            dplyr::mutate(RecoveryPercentage = round(Recovery * 100, 0)) |>  # Round Recovery to 0 decimals
            dplyr::select(`Replicate Name`, `Sample Type`, RecoveryPercentage)  # Select only relevant columns

        # Replace NA values with 0 if necessary
        RECOVERY[is.na(RECOVERY)] <- 0

        return(RECOVERY)
    })

    # Render the QA/QC datatable with only RecoveryPercentage column
    output$table2 <- DT::renderDT({
        RECOVERY <- Skyline_recovery()  # Get the reactive data

        # Check if RECOVERY has data
        req(nrow(RECOVERY) > 0)

        DT::datatable(RECOVERY,
                      filter = "top",
                      extensions = c("Buttons", "Scroller"),
                      options = list(
                          scrollY = 650,
                          scrollX = 500,
                          deferRender = TRUE,
                          scroller = TRUE,
                          buttons = list(
                              list(
                                  extend = "excel",
                                  filename = "Samples_recovery",
                                  title = NULL,
                                  exportOptions = list(modifier = list(page = "all"))
                              ),
                              list(
                                  extend = "csv",
                                  filename = "Samples_recovery",
                                  title = NULL,
                                  exportOptions = list(modifier = list(page = "all"))
                              ),
                              list(
                                  extend = "colvis",
                                  targets = 0,
                                  visible = FALSE
                              )
                          ),
                          dom = "lBfrtip",
                          fixedColumns = TRUE
                      ),
                      rownames = FALSE
        )
    })


    ###################################################################LOD####################################

    # Define a reactive block for the LOD table
    LOD_summary <- reactive({
        # Ensure data is available
        df_samples <- ConcentrationB  # Use the reactive data source for Samples_Concentration

        req(df_samples)

        # Filter for 'Blank' Sample Type and calculate the standard deviation of 'Samples_Concentration'
        lod_sd <- df_samples |>
            dplyr::filter('Sample Type' == "Blank") |>
            dplyr::summarize(LOD = sd(Concentration, na.rm = TRUE) * 3)  # Multiply by 3 for LOD

        # Convert to data frame
        lod_sd <- as.data.frame(lod_sd)

        lod_sd
    })

    # Render the LOD table
    output$LOD <- DT::renderDT({
        lod_data <- LOD_summary()  # Get the reactive data

        DT::datatable(lod_data,
                      options = list(
                          pageLength = 10,
                          dom = 't'  # Table only, without additional controls
                      ),
                      rownames = FALSE
        )
    })


    # Close the app when the session ends
    if(!interactive()) {
        session$onSessionEnded(function() {
            stopApp()
            q("no")
        })
    }

}



# Run the application
shinyApp(ui = ui, server = server)

