#' CPquant: Shiny CPs Quantification for Skyline Output
#' @param ...
#'
#' @import readxl
#' @import nnls
#' @import crosstalk
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import shiny
#' @import DT
#' @import plotly
#' @import purrr
#' @import markdown
#' @importFrom stats lm


CPquant <- function(...){

    options(shiny.maxRequestSize = 500 * 1024^2)

    ui <- shiny::navbarPage("Quantification by deconvolution from Skyline output",
                            shiny::tabPanel("Quantification Inputs",
                                            shiny::fluidPage(shiny::sidebarLayout(
                                                shiny::sidebarPanel(
                                                    width = 3,
                                                    shiny::fileInput("fileInput", "Import excel file from Skyline",
                                                                     accept = c('xlsx')),
                                                    shiny::radioButtons("blankSubtraction",
                                                                        label = "Subtraction with blank?",
                                                                        choices = c("Yes, by avg area of blanks", "No"), selected = "No"),
                                                    shiny::radioButtons("correctWithRS", label = "Correct with RS?",
                                                                        choices = c("Yes", "No"), selected = "No"),
                                                    shiny::radioButtons("calculateRecovery",
                                                                        label = "Calculate recovery? (req QC samples)",
                                                                        choices = c("Yes", "No"), selected = "No"),
                                                    shiny::radioButtons("calculateMDL",
                                                                        label = "Calculate MDL? (req blank samples)",
                                                                        choices = c("Yes", "No"), selected = "No"),
                                                    shiny::radioButtons("standardTypes", label = "Types of standards",
                                                                        choices = c("Mixtures"), selected = "Mixtures"), #will add single chain std later
                                                    shiny::actionButton('go', 'Proceed', width = "100%"),
                                                    shiny::uiOutput("defineVariables")
                                                ),
                                                shiny::mainPanel(
                                                    width = 9,
                                                    plotly::plotlyOutput("plot_Skyline_output"),
                                                    DT::DTOutput("table_Skyline_output")

                                                )
                                            )
                                            )),
                            shiny::tabPanel(
                                "Input summary",
                                shiny::sidebarPanel(
                                    width = 2, # max 12
                                    shiny::radioButtons("navSummary", "Choose tab:",
                                                        choices = c("Std Calibration Curves", "Quan to Qual ratio"),
                                                        selected = "Std Calibration Curves")
                                ),
                                shiny::mainPanel(
                                    width = 10,
                                    shiny::conditionalPanel(
                                        condition = "input.navSummary == 'Std Calibration Curves'",
                                        plotly::plotlyOutput("CalibrationSCCPs"),
                                        shiny::tags$br(),
                                        plotly::plotlyOutput("CalibrationMCCPs"),
                                        shiny::tags$br(),
                                        plotly::plotlyOutput("CalibrationLCCPs")
                                    ),
                                    shiny::conditionalPanel(
                                        condition = "input.navSummary == 'Quan to Qual ratio'",
                                        plotly::plotlyOutput("RatioQuantToQual")
                                    )
                                )
                            ),

                            shiny::tabPanel(
                                "Quantification summary",
                                shiny::fluidPage(DT::DTOutput("quantTable"))
                            ),
                            shiny::tabPanel(
                                "QA/QC",
                                shiny::mainPanel(
                                    DT::DTOutput("table_recovery"),   # First table output (Skyline recovery data)
                                    br(),                     # Optional line break to add space between the tables
                                    DT::DTOutput("table_MDL")       # Second table output (LOD table)
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

        Skyline_output <- reactive({
            req(input$fileInput) #requires that the input is available

            # Create a Progress object
            progress <- shiny::Progress$new()
            on.exit(progress$close())

            progress$set(message = "WAIT! Loading data...", value = 0)

            # Read the Excel file
            progress$set(value = 0.3, detail = "Reading Excel file")
            df <- readxl::read_excel(input$fileInput$datapath)

            progress$set(value = 0.6, detail = "Processing data")
            df <- df |>
                dplyr::mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) |>
                dplyr::mutate(Area = as.numeric(Area)) |>
                dplyr::mutate(Area = replace_na(Area, 0)) |>
                dplyr::mutate(C_part = stringr::str_extract(Molecule, "C\\d+"),
                              Cl_part = stringr::str_extract(Molecule, "Cl\\d+"),
                              C_number = as.numeric(stringr::str_extract(C_part, "\\d+")),
                              Cl_number = as.numeric(stringr::str_extract(Cl_part, "\\d+")),
                              PCA = stringr::str_c(C_part, Cl_part, sep = ""))

            progress$set(value = 0.8, detail = "Applying corrections")

            #  Normalize data based on 'Correct with RS' input
            if (input$correctWithRS == "Yes" & any(df$`Molecule List` == "RS")){
                df <- df |>
                    dplyr::group_by(`Replicate Name`) |>
                    dplyr::mutate(Area = Area / first(Area[`Molecule List`== "RS" & `Isotope Label Type` == "Quan"])) |>
                    dplyr::ungroup()
            }

            # Calculate the average blank value
            if (input$blankSubtraction == "Yes, by avg area of blanks"){
                df_blank <- df |>
                    dplyr::filter(`Sample Type` == "Blank") |>
                    dplyr::group_by(Molecule, `Molecule List`, `Isotope Label Type`) |>
                    dplyr::summarize(AverageBlank = mean(Area, na.rm = TRUE)) |>
                    dplyr::ungroup() |>
                    dplyr::filter(!`Molecule List` %in% c("IS", "RS", "VS"))

                df <- df |>
                    dplyr::full_join(df_blank) |>
                    dplyr::mutate(AverageBlank = tidyr::replace_na(AverageBlank, 0)) |>
                    dplyr::mutate(Area = dplyr::case_when(`Sample Type` == "Unknown" ~ Area - AverageBlank, .default = Area)) |>
                    dplyr::mutate(Area = ifelse(Area <0, 0, Area))
            }

            progress$set(value = 1, detail = "Complete")

            return(df)
        })

        # defineVariablesUI in separate file UI_components.R
        output$defineVariables <- shiny::renderUI({
            defineVariablesUI(Skyline_output())
        })


        # Set reactive values from user input

        removeRsquared <- shiny::eventReactive(input$go, {as.numeric(input$removeRsquared)})
        standardAnnoColumn <- shiny::eventReactive(input$go, {as.character(input$standardAnnoColumn)})
        removeSamples <- shiny::eventReactive(input$go, {as.character(input$removeSamples)})


        #Render raw table
        output$table_Skyline_output <- DT::renderDT({
            DT::datatable(Skyline_output(),
                          options = list(
                              paging = TRUE,
                              pageLength = 50
                          )
            )
        })


        # Render reactive summary statistics and plots of raw input BEFORE quantification


        output$plot_Skyline_output <- plotly::renderPlotly({
            Skyline_output() |>
                dplyr::filter(`Isotope Label Type` == "Quan") |>
                dplyr::mutate(OrderedMolecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)]))) |>  # Create a composite ordering factor

                plotly::plot_ly(
                    x = ~ OrderedMolecule,
                    y = ~ Area,
                    color = ~ `Sample Type`,
                    type = "box",
                    text = ~paste(
                        "Homologue: ", PCA,
                        "<br>Sample: ", `Replicate Name`,
                        "<br>Area:", round(Area, 2)
                    ),
                    hoverinfo = "text"
                )
        })



        #----------------START: Deconvolution script------------------#

        shiny::observeEvent(input$go, {

            progress <- shiny::Progress$new()
            on.exit(progress$close())

            progress$set(message = "WAIT! Processing data...", value = 0)

            # remove samples if selected by removeSamples input

            if(!is.null(removeSamples()) && length(removeSamples()) > 0){
                Skyline_output_filt <- Skyline_output() |>
                    dplyr::filter(!`Replicate Name` %in% removeSamples())
            } else{
                Skyline_output_filt <- Skyline_output()
            }



            ##### PREPARE FOR DECONVOLUTION #######
            # Prepare for deconvolution for standards
            progress$set(value = 0.2, detail = "Preparing standards data")
            CPs_standards <- Skyline_output_filt |>
                dplyr::filter(`Sample Type` == "Standard", #make sure the stds are not blank corrected
                              !`Molecule List` %in% c("IS", "RS", "VS"), # dont include IS, RS, VS
                              `Isotope Label Type` == "Quan", # use only Quant ions
                              !!dplyr::sym(standardAnnoColumn()) != "NA") |>
                #dplyr::group_by(!!dplyr::sym(standardAnnoColumn()), Molecule, `Molecule List`, C_number, Cl_number, PCA) |> #"Batch Name" is default. !!sym(standardAnnoColumn) unquotes the string variable and converts it to a symbol that dplyr can understand within the group_by() function
                dplyr::group_by(!!dplyr::sym(standardAnnoColumn()), `Sample Type`, Molecule, `Molecule List`, C_number, Cl_number, PCA) |> #"Batch Name" is default. !!sym(standardAnnoColumn) unquotes the string variable and converts it to a symbol that dplyr can understand within the group_by() function
                tidyr::nest() |>
                dplyr::mutate(models = purrr::map(data, ~lm(Area ~ `Analyte Concentration`, data = .x))) |>
                dplyr::mutate(coef = purrr::map(models, coef)) |>
                dplyr::mutate(Response_factor = purrr::map_dbl(models, ~ coef(.x)["`Analyte Concentration`"]))|> #get the slope which will be RF
                dplyr::mutate(intercept = purrr::map(coef, purrr::pluck("(Intercept)"))) |>
                dplyr::mutate(rsquared = purrr::map(models, summary)) |> #first create a data frame list with the model
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







            ###### This is for mixtures, single chain stds will be added later

            progress$set(value = 0.4, detail = "Processing standards")
            if(input$standardTypes == "Mixtures"){
                # For SCCPs
                CPs_standards_S <- CPs_standards  |>
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



                ##### ADD from standards calibration (plots.R)
                output$CalibrationSCCPs <- plotly::renderPlotly({
                    plot_cal_SCCPs(CPs_standards_S, standardAnnoColumn())
                })

                output$CalibrationMCCPs <- plotly::renderPlotly({
                    plot_cal_MCCPs(CPs_standards_M, standardAnnoColumn())
                })

                output$CalibrationLCCPs <- plotly::renderPlotly({
                    plot_cal_LCCPs(CPs_standards_L, standardAnnoColumn())
                })

            }






            # Prepare for deconvolution of samples


            progress$set(value = 0.6, detail = "Preparing sample data")
            CPs_samples <- Skyline_output_filt |>
                dplyr::filter(
                    #`Sample Type` == "Unknown",
                    `Sample Type` %in% c("Unknown", "Blank"), #include both unknown and blank
                    !`Molecule List` %in% c("IS", "RS", "VS"), # remove IS, RS, VS
                    `Isotope Label Type` == "Quan") |>
                dplyr::group_by(`Replicate Name`) |>  # Group by Replicate Nam
                dplyr::mutate(Relative_distribution = Area / sum(Area, na.rm = TRUE)) |>
                dplyr::ungroup() |>
                #dplyr::select(`Replicate Name`, `Sample Type`, Molecule, Area, Relative_distribution) |>
                dplyr::select(-`Mass Error PPM`, -`Isotope Label Type`, -`Chromatogram Precursor M/Z`, -`Analyte Concentration`, -`Batch Name`) |>
                dplyr::mutate(dplyr::across(Relative_distribution, ~replace(., is.nan(.), 0)))  # Replace NaN with zero




            ############################################################################### DECONVOLUTION #############################################################################

            progress$set(value = 0.8, detail = "Performing deconvolution")
            # Deconvolution process
            CPs_standards_input <- CPs_standards |>
                dplyr::select(Molecule, !!dplyr::sym(standardAnnoColumn()), Response_factor) |>
                tidyr::pivot_wider(names_from = !!dplyr::sym(standardAnnoColumn()), values_from = "Response_factor")


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


            # Ensure combined_standard is correctly defined as a matrix prior to the deconvolution

            # First populate combined_standard with CPs_standards_input
            combined_standard <- CPs_standards_input  |>
                tibble::column_to_rownames(var = "Molecule") |>
                #dplyr::select(-Molecule) |>
                as.matrix()


            # Ensure combined_sample is correctly defined with nested data frames prior to the deconvolution
            combined_sample <- CPs_samples  |>
                dplyr::group_by(`Replicate Name`, `Sample Type`) |>
                tidyr::nest() |>
                dplyr::ungroup()



            # Apply the perform_deconvolution function to each nested data frame
            deconvolution <- combined_sample |>
                #perform_deconvolution on only Relative_distribution in the nested data frame
                dplyr::mutate(result = purrr::map(data, ~ perform_deconvolution(dplyr::select(.x, Relative_distribution), combined_standard))) |>
                dplyr::mutate(total_sum = purrr::map_dbl(result, ~sum(.x$deconv_resolved))) |>
                dplyr::mutate(deconv_coef = purrr::map(result, ~as_tibble(list(deconv_coef = .x$deconv_coef, `Batch Name` = names(.x$deconv_coef))))) |>
                dplyr::mutate(rsquared = purrr::map(result, purrr::pluck("r_squared"))) |>
                dplyr::mutate(deconv_resolved = purrr::map(result, ~tibble::as_tibble(list(deconv_resolved = .x$deconv_resolved, Molecule = rownames(.x$deconv_resolved))))) |>
                #calculate the Concentration by multiplying the total_sum with Relative_distribution inside data nested object
                dplyr::mutate(data = purrr::map2(data, total_sum, ~dplyr::mutate(.x, Concentration = .y *Relative_distribution))) |>
                dplyr::select(-result)



            ########################################################## Calculate the concentration in ng/uL ###############################################################


            progress$set(value = 0.9, detail = "Calculating final results")

                Samples_Concentration <- deconvolution |>
                dplyr::select(`Replicate Name`, `Sample Type`, data) |>
                dplyr::mutate(data = purrr::map(data, ~ select(.x, Molecule, Concentration))) |>
                unnest(data) |>
                pivot_wider(names_from = Molecule, values_from = Concentration)




            # Extract deconv_coef from results and create a new data frame
            # deconv_coef_df <- Deconvolution |>
            #     dplyr::mutate(deconv_coef = purrr::map(result, "deconv_coef")) |>
            #     dplyr::mutate(rsquared = purrr::map(result, purrr::pluck("r_squared"))) |> # then pluck only the r.squared value
            #     dplyr::select(`Replicate Name`, `Sample Type`, rsquared,  deconv_coef) |>
            #     tidyr::unnest_wider(deconv_coef, names_sep = "_") |>
            #     dplyr::rename_with(~ str_replace_all(., "deconv_coef_", ""))
            #
            # #Calculate the response of the standards
            #
            # #Remove the replicate name to generate vectors:
            # deconv_coef_df_matrix<- deconv_coef_df |>
            #     #tibble::column_to_rownames(var = "Replicate Name")
            #     select(-`Replicate Name`, -`Sample Type`, -rsquared)
            #
            #
            #
            # # Initialize an empty list to store results
            # result_list <- list()
            #
            # # Iterate through each row of deconv_coef_df_matrix
            # for (i in 1:nrow(deconv_coef_df_matrix)) {
            #
            #     # Extract row vector from deconv_coef_df_matrix
            #     deconv_coef_vector <- as.numeric(deconv_coef_df_matrix[i, ])
            #
            #     # # I had this one before, but to make sure
            #     # combined_standard <- CPs_standards_input |>
            #     #     select(-Molecule)  |>
            #     #     mutate(across(everything(), as.numeric)) |>
            #     #     as.matrix()
            #
            #     # Perform element-wise multiplication
            #     result <- sweep(combined_standard, 2, deconv_coef_vector, `*`)
            #
            #     # Create data frame with column names from CPs_standards_input
            #     result_df <- as.data.frame(result)
            #     colnames(result_df) <- colnames(CPs_standards_input)[-which(names(CPs_standards_input) == "Molecule")]
            #
            #     # Assign name to the result_df from deconv_coef_df
            #     replicate_name <- deconv_coef_df$`Replicate Name`[i]
            #
            #     # Store result in result_list with the corresponding name
            #     result_list[[replicate_name]] <- result_df
            # }
            #
            #
            # # Combine all data frames into a single data frame with 'Replicate Name'
            # final_df <- do.call(rbind, Map(function(df, name) {
            #     df$`Replicate Name` <- name
            #     df <- df[, c("Replicate Name", setdiff(names(df), "Replicate Name"))]
            #     df
            # }, result_list, names(result_list)))
            #
            # # Add CPs_standards_input$Molecule column to final_df
            # final_df$Molecule <- CPs_standards_input$Molecule
            #
            #
            #
            # #Organize the data
            # final_df_tidy <- final_df|>
            #     group_by(`Replicate Name`) |>
            #     nest() |>
            #     ungroup()
            #
            #
            # #Total sum the values for each replicate
            #
            # #Remove the molecule
            # final_df_matrix<- final_df |>
            #     select(-Molecule) |>
            #     group_by(`Replicate Name`) |>
            #     nest()
            #
            # # Initialize an empty data frame to store results
            # total_sums_df <- data.frame(
            #     `Replicate Name` = character(),
            #     `Total Sum` = numeric(),
            #     stringsAsFactors = FALSE)
            #
            # # Iterate through each row of final_df_grouped
            # for (i in 1:nrow(final_df_matrix)) {
            #     # Extract nested data frame
            #     nested_df <- final_df_matrix$data[[i]]
            #
            #     # Calculate total sum for the current `Replicate Name`
            #     `Replicate Name` <- final_df_matrix$`Replicate Name`[[i]]
            #     total_sum <- sum(colSums(nested_df[, -1]))  # Exclude the grouping column
            #
            #     # Append results to total_sums_df
            #     total_sums_df <- rbind(total_sums_df, data.frame(
            #         `Replicate Name` = `Replicate Name`,
            #         `Total Sum` = total_sum
            #     ))
            # }
            #
            #
            #
            #
            #
            #
            # ################################################### FINAL RESULTS ####################################################################
            #
            #
            # # # Perform operations to reorganize data
            # # Final_results <- merged_df |>
            # #     dplyr::select(-Area, -Relative_distribution, -Calculated_RF, -Measured_Area, -Concentration) |> # Remove unwanted columns
            # #     tidyr::pivot_wider(
            # #         names_from = Replicate.Name,  # Use values from Replicate.Name as new column names
            # #         values_from = ConcentrationDetailed  # Use values from ConcentrationDetailed to fill the new columns
            # #     )
            #
            # # Merge total_sums_df into CPs_samples based on Replicate Name
            #
            # total_sums_df <- total_sums_df |>  mutate(`Replicate Name` = Replicate.Name)
            #
            # Concentration <- CPs_samples  |>
            #     left_join(total_sums_df, by = "Replicate Name")  |>
            #     mutate(Concentration = `Relative_distribution` * `Total.Sum`)
            #
            #
            # Concentration<-Concentration |>
            #     group_by(`Replicate Name`) |>
            #     distinct( `Molecule`, Concentration) |>
            #     nest()
            #
            # # Perform operations to reorganize data
            # reorganized_data <- Concentration  |>
            #     unnest(cols = c(data)) |>
            #     distinct(`Replicate Name`, Molecule, .keep_all = TRUE)  |>
            #     pivot_wider(names_from = Molecule, values_from = Concentration)
            # reorganized_data <- t(reorganized_data) #transpose
            #
            # #Make the first row (replicate names) the column names
            # colnames(reorganized_data) <- reorganized_data[1, ]
            # Samples_Concentration <- reorganized_data[-1, ]
            # # Convert the result back to a data frame
            # Samples_Concentration <- as.data.frame(Samples_Concentration)
            # Samples_Concentration<- Samples_Concentration |>
            #     mutate(Molecule = CPs_samples_input$Molecule)|>
            #     relocate(Molecule, .before = everything())



            progress$set(value = 1, detail = "Complete")

            ### END: Deconvolution script


            # Render table
            output$quantTable <- DT::renderDT({
                DT::datatable(
                    #Final_results,
                    Samples_Concentration,
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



            ###########################################################RECOVERY#######################################################

            if(input$calculateRecovery == "Yes") {
                # Recovery calculations
                recovery_data <- Skyline_output_filt |>
                    dplyr::filter(`Isotope Label Type` == "Quan",
                                  `Molecule List` %in% c("IS", "RS"))


                RECOVERY <- recovery_data |>  # Calculate recovery
                    tidyr::pivot_wider(
                        id_cols = c(`Replicate Name`, `Sample Type`),
                        names_from = `Molecule List`,
                        values_from = Area
                    ) |>
                    dplyr::mutate(across(c(IS, RS), ~replace_na(.x, 0)))

                # Calculate QC ratio
                qc_ratio <- RECOVERY |>
                    dplyr::filter(`Sample Type` == "Quality Control") |>
                    dplyr::mutate(RatioStd = IS / RS) |>
                    dplyr::summarize(AverageRatio = mean(RatioStd, na.rm = TRUE))

                # Calculate sample recovery
                RECOVERY <- RECOVERY |>
                    dplyr::filter(`Sample Type` %in% c("Unknown", "Blank")) |>
                    dplyr::mutate(
                        RatioSample = IS / RS,
                        Recovery = RatioSample / as.numeric(qc_ratio$AverageRatio),
                        RecoveryPercentage = round(Recovery * 100, 0)
                    ) |>
                    dplyr::select(`Replicate Name`, `Sample Type`, RecoveryPercentage)


                # Render recovery table
                output$table_recovery <- DT::renderDT({
                    DT::datatable(RECOVERY,
                                  filter = "top",
                                  extensions = c("Buttons", "Scroller"),
                                  options = list(
                                      scrollY = 650,
                                      scrollX = 500,
                                      deferRender = TRUE,
                                      scroller = TRUE,
                                      buttons = list(
                                          list(extend = "excel",
                                               filename = "Samples_recovery",
                                               title = NULL,
                                               exportOptions = list(modifier = list(page = "all"))),
                                          list(extend = "csv",
                                               filename = "Samples_recovery",
                                               title = NULL,
                                               exportOptions = list(modifier = list(page = "all"))),
                                          list(extend = "colvis",
                                               targets = 0,
                                               visible = FALSE)
                                      ),
                                      dom = "lBfrtip",
                                      fixedColumns = TRUE
                                  ),
                                  rownames = FALSE)
                })
            }


            if(input$calculateMDL == "Yes") {
            #MDL calculations (need to take into account if blank subtraction affect or not)
            MDL_data <- Samples_Concentration |>
                dplyr::filter(`Sample Type` == "Blank") |>
                dplyr::rowwise() %>%
                dplyr::mutate(sum_PCA = sum(dplyr::c_across(-c(`Replicate Name`, `Sample Type`)), na.rm = TRUE)) %>%
                dplyr::ungroup() |>
                dplyr::summarize(
                    MDL_sumPCA = 3 * sd(sum_PCA, na.rm = TRUE),
                    number_of_blanks = n()
                )
            }



            # Render MDL table
            output$table_MDL <- DT::renderDT({
                DT::datatable(MDL_data,
                              options = list(
                                  pageLength = 10,
                                  dom = 't'
                              ),
                              rownames = FALSE)
            })

browser()

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

}
