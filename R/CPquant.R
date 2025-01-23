#' CPquant: Shiny CPs Quantification for Skyline Output
#' @param ...
#'
#' @import shiny
#' @import htmlwidgets
#' @import ggplot2
#' @import readxl
#' @import nnls
#' @import crosstalk
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import DT
#' @import plotly
#' @import purrr
#' @import markdown
#' @import openxlsx
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
                                                    shiny::textInput("quantificationUnit", "Enter Quantification unit:"),
                                                    shiny::radioButtons("blankSubtraction",
                                                                        label = "Subtraction with blank?",
                                                                        choices = c("Yes, by avg area of blanks", "No"), selected = "No"),
                                                    shiny::radioButtons("correctWithRS", label = "Correct with RS area?",
                                                                        choices = c("Yes", "No"), selected = "No"),
                                                    shiny::radioButtons("calculateRecovery",
                                                                        label = "Calculate recovery? (req QC samples)",
                                                                        choices = c("Yes", "No"), selected = "No"),
                                                    shiny::radioButtons("calculateMDL",
                                                                        label = "Calculate MDL? (req blank samples)",
                                                                        choices = c("Yes", "No"), selected = "No"),
                                                    shiny::radioButtons("standardTypes", label = "Types of standards",
                                                                        choices = c("Group Mixtures"), selected = "Group Mixtures"), #will add single chain std later
                                                    shiny::tags$div(
                                                        title = "Wait until import file is fully loaded before pressing!",
                                                        shiny::actionButton('go', 'Proceed', width = "100%")
                                                    ),
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
                                        plotly::plotlyOutput("CalibrationCurves")
                                    ),
                                    shiny::conditionalPanel(
                                        condition = "input.navSummary == 'Quan to Qual ratio'",
                                        plotly::plotlyOutput("RatioQuantToQual", height = "80vh", width = "100%")
                                    )
                                )
                            ),

                            shiny::tabPanel(
                                "Quantification summary",
                                shiny::fluidPage(
                                    downloadButton("downloadResults", "Export all results to Excel"),
                                    shiny::tags$br(),
                                    DT::DTOutput("quantTable"),
                                    shiny::tags$br(),
                                    plotly::plotlyOutput("sampleContributionPlot")
                                )
                            ),
                            shiny::tabPanel(
                                "Homologue Group Patterns",
                                fluidRow(
                                    column(
                                        width = 2,
                                        shiny::radioButtons("plotHomologueGroups", "Choose tab:",
                                                            choices = c("All Samples Overview", "Samples Overlay", "Samples Panels"),
                                                            selected = "All Samples Overview"),
                                        shiny::conditionalPanel(
                                            condition = "input.plotHomologueGroups == 'Samples Overlay'",
                                            shiny::uiOutput("sampleSelectionUIOverlay")
                                        ),
                                        shiny::conditionalPanel(
                                            condition = "input.plotHomologueGroups == 'Samples Panels'",
                                            shiny::uiOutput("sampleSelectionUIComparisons")
                                        ),
                                        shiny::tags$div(
                                            title = "WAIT after pressing..Might take some time before plot shows!",
                                            shiny::actionButton('go2', 'Plot', width = "100%")
                                        )
                                    ),
                                    column(
                                        width = 10,
                                        shiny::conditionalPanel(
                                            condition = "input.plotHomologueGroups == 'All Samples Overview'",
                                            shiny::plotOutput("plotHomologuePatternStatic", height = "80vh", width = "100%")
                                        ),
                                        shiny::conditionalPanel(
                                            condition = "input.plotHomologueGroups == 'Samples Overlay'",
                                            plotly::plotlyOutput("plotHomologuePatternOverlay", height = "80vh", width = "100%")
                                        ),
                                        shiny::conditionalPanel(
                                            condition = "input.plotHomologueGroups == 'Samples Panels'",
                                            plotly::plotlyOutput("plotHomologuePatternComparisons", height = "80vh", width = "100%")
                                        )
                                    )
                                )
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
            df <- readxl::read_excel(input$fileInput$datapath) #outputs a tibble

            progress$set(value = 0.6, detail = "Processing data")
            df <- df |>
                dplyr::rename(Replicate_Name = `Replicate Name`) |>
                dplyr::rename(Sample_Type = `Sample Type`) |>
                dplyr::rename(Molecule_List = `Molecule List`) |>
                dplyr::rename(Mass_Error_PPM = `Mass Error PPM`) |>
                dplyr::rename(Isotope_Label_Type = `Isotope Label Type`) |>
                dplyr::rename(Chromatogram_Precursor_MZ = `Chromatogram Precursor M/Z`) |>
                dplyr::rename(Analyte_Concentration = `Analyte Concentration`) |>
                dplyr::rename(Batch_Name = `Batch Name`) |>
                dplyr::mutate(Analyte_Concentration = as.numeric(Analyte_Concentration)) |>
                dplyr::mutate(Area = as.numeric(Area)) |>
                dplyr::mutate(Area = replace_na(Area, 0)) |>
                dplyr::mutate(C_homologue = stringr::str_extract(Molecule, "C\\d+"),
                              Cl_homologue = stringr::str_extract(Molecule, "Cl\\d+"),
                              C_number = as.numeric(stringr::str_extract(C_homologue, "\\d+")),
                              Cl_number = as.numeric(stringr::str_extract(Cl_homologue, "\\d+")),
                              PCA = stringr::str_c(C_homologue, Cl_homologue, sep = ""))

            progress$set(value = 0.8, detail = "Applying corrections")

            #  Normalize data based on 'Correct with RS' input
            if (input$correctWithRS == "Yes" & any(df$Molecule_List == "RS")){
                df <- df |>
                    dplyr::group_by(Replicate_Name) |>
                    dplyr::mutate(Area = Area / first(Area[Molecule_List== "RS" & Isotope_Label_Type == "Quan"])) |>
                    dplyr::ungroup()
            }

            # Calculate the average blank value
            if (input$blankSubtraction == "Yes, by avg area of blanks"){
                df_blank <- df |>
                    dplyr::filter(Sample_Type == "Blank") |>
                    dplyr::group_by(Molecule, Molecule_List, Isotope_Label_Type) |>
                    dplyr::summarize(AverageBlank = mean(Area, na.rm = TRUE)) |>
                    dplyr::ungroup() |>
                    dplyr::filter(!Molecule_List %in% c("IS", "RS", "VS"))

                df <- df |>
                    dplyr::full_join(df_blank) |>
                    dplyr::mutate(AverageBlank = tidyr::replace_na(AverageBlank, 0)) |>
                    dplyr::mutate(Area = dplyr::case_when(Sample_Type == "Unknown" ~ Area - AverageBlank, .default = Area)) |> #only blank subtraction of Unknown Sample Type
                    dplyr::mutate(Area = ifelse(Area <0, 0, Area)) #replace negative Area with 0 after blank subtraction
            }


            if (input$standardTypes == "Group Mixtures") {
                df <- df |>
                    dplyr::mutate(Quantification_Group = stringr::str_extract(Batch_Name, "^[^_]+"))
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
        #standardAnnoColumn <- shiny::eventReactive(input$go, {as.character(input$standardAnnoColumn)})
        removeSamples <- shiny::eventReactive(input$go, {as.character(input$removeSamples)})
        Samples_Concentration <- reactiveVal() # Create a reactive value to store Samples_Concentration after deconvolution


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
                dplyr::filter(Isotope_Label_Type == "Quan") |>
                dplyr::mutate(OrderedMolecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)]))) |>  # Create a composite ordering factor

                plotly::plot_ly(
                    x = ~ OrderedMolecule,
                    y = ~ Area,
                    color = ~ Sample_Type,
                    type = "box",
                    text = ~paste(
                        "Homologue: ", PCA,
                        "<br>Sample: ", Replicate_Name,
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
                    dplyr::filter(!Replicate_Name %in% removeSamples())
            } else{
                Skyline_output_filt <- Skyline_output()
            }



            ##### PREPARE FOR DECONVOLUTION #######
            # Prepare for deconvolution for standards
            progress$set(value = 0.2, detail = "Preparing standards data")

            # Prepare data frame for all standards to be used in deconvolution
            if(input$standardTypes == "Group Mixtures"){

                CPs_standards <- Skyline_output_filt |>
                    dplyr::filter(Sample_Type == "Standard",
                                  !Molecule_List %in% c("IS", "RS", "VS"), # dont include IS, RS, VS
                                  Isotope_Label_Type == "Quan", # use only Quant ions
                                  Batch_Name != "NA") |>
                    dplyr::mutate(
                        C_range = stringr::str_extract_all(Quantification_Group, "\\d+"),
                        C_min = as.numeric(purrr::map_chr(C_range, ~.x[1])),
                        C_max = as.numeric(purrr::map_chr(C_range, function(x) {if(length(x) > 1) {x[2]} else {x[1]}}))) |>
                    dplyr::select(-C_range) |>
                    dplyr::group_by(Batch_Name, Sample_Type, Molecule, Molecule_List, C_number, Cl_number, PCA, Quantification_Group, C_min, C_max) |> #Batch_Name is default. !!sym(standardAnnoColumn) unquotes the string variable and converts it to a symbol that dplyr can understand within the group_by() function
                    tidyr::nest() |>
                    dplyr::filter(C_number >= C_min & C_number <= C_max) |>
                    dplyr::mutate(models = purrr::map(data, ~lm(Area ~ Analyte_Concentration, data = .x))) |>
                    dplyr::mutate(coef = purrr::map(models, coef)) |>
                    dplyr::mutate(RF = purrr::map_dbl(models, ~ coef(.x)["Analyte_Concentration"]))|> #get the slope which will be the RF
                    dplyr::mutate(intercept = purrr::map(coef, purrr::pluck("(Intercept)"))) |>
                    dplyr::mutate(rsquared = purrr::map(models, summary)) |> #first create a data frame list with the model
                    dplyr::mutate(rsquared = purrr::map(rsquared, purrr::pluck("r.squared"))) |> # then pluck only the r.squared value
                    dplyr::select(-coef) |>  # remove coef variable since it has already been plucked
                    tidyr::unnest(c(RF, intercept, rsquared)) |>  #removing the list type for these variables
                    dplyr::mutate(RF = if_else(RF < 0, 0, RF)) |> # replace negative RF with 0
                    dplyr::mutate(rsquared = ifelse(is.nan(rsquared), 0, rsquared)) |>
                    dplyr::mutate(RF = if_else(rsquared < removeRsquared(), 0, RF)) |> #keep RF only if rsquared is above removeRsquared input
                    dplyr::ungroup() |>
                    dplyr::group_by(Batch_Name) |> #grouping by the standards
                    dplyr::mutate(Sum_RF_group = sum(RF, na.rm = TRUE)) |>
                    dplyr::ungroup()




                # Prepare for deconvolution of samples
                progress$set(value = 0.6, detail = "Preparing sample data")

                # Character vector for unique Batch Name
                #unique_chars <- unique(stringr::str_extract(CPs_standards$Batch_Name, "^[^_]+"))

                CPs_samples <- Skyline_output_filt |>
                    dplyr::filter(
                        #Sample_Type == "Unknown",
                        Sample_Type %in% c("Unknown", "Blank"), #include both unknown and blank
                        !Molecule_List %in% c("IS", "RS", "VS"), # remove IS, RS, VS
                        Isotope_Label_Type == "Quan") |>
                    # Quantification_Group
                    dplyr::group_by(Replicate_Name) |>  # Group by Replicate Name
                    dplyr::mutate(Relative_Area = Area / sum(Area, na.rm = TRUE)) |> #Relative area
                    dplyr::ungroup() |>
                    dplyr::select(-Mass_Error_PPM, -Isotope_Label_Type, -Chromatogram_Precursor_MZ, -Analyte_Concentration, -Batch_Name) |>
                    dplyr::mutate(dplyr::across(Relative_Area, ~replace(., is.nan(.), 0)))  # Replace NaN with zero


                CPs_samples_input <- CPs_samples |>
                    dplyr::select(Molecule, Replicate_Name, Relative_Area) |>
                    tidyr::pivot_wider(names_from = "Replicate_Name", values_from = "Relative_Area")



                # Ensure combined_sample is correctly defined with nested data frames prior to deconvolution
                combined_sample <- CPs_samples  |>
                    dplyr::group_by(Replicate_Name, Sample_Type) |>
                    tidyr::nest() |>
                    dplyr::ungroup()


                ###### Plot the calibration curves ######

                output$CalibrationCurves <- plotly::renderPlotly({
                    plot_calibration_curves(CPs_standards)
                })

                output$RatioQuantToQual <- plotly::renderPlotly({

                    Skyline_output_filt |>
                        dplyr::group_by(Replicate_Name, Molecule) |>
                        dplyr::mutate(Quan_Area = ifelse(Isotope_Label_Type == "Quan", Area, NA)) |>
                        tidyr::fill(Quan_Area, .direction = "downup") |>
                        dplyr::mutate(QuanMZ = ifelse(Isotope_Label_Type == "Quan", Chromatogram_Precursor_MZ, NA)) |>
                        tidyr::fill(QuanMZ, .direction = "downup") |>
                        dplyr::mutate(QuanQualRatio = ifelse(Isotope_Label_Type == "Qual", Quan_Area/Area, 1)) |>
                        tidyr::replace_na(list(QuanQualRatio = 0)) |>
                        dplyr::mutate(QuanQualMZ = paste0(QuanMZ,"/",Chromatogram_Precursor_MZ)) |>
                        dplyr::ungroup() |>
                        dplyr::select(Replicate_Name, Sample_Type, Molecule_List, Molecule, QuanQualMZ, QuanQualRatio) |>
                        plot_ly(x = ~Replicate_Name, y = ~QuanQualRatio, type = 'violin', color = ~Sample_Type,
                                text = ~paste(Replicate_Name, "<br>", Molecule_List, "<br>", Molecule,
                                              "<br>Quan/Qual MZ: ", QuanQualMZ,
                                              "<br>Ratio: ", round(QuanQualRatio, 2)),
                                hoverinfo = "text") |>
                        layout(title = 'Quan-to-Qual Ratio',
                               xaxis = list(title = 'Replicate Name'),
                               yaxis = list(title = 'Quan-to-Qual Ratio'))
                })


                #### DECONVOLUTION ####

                progress$set(value = 0.8, detail = "Performing deconvolution")


                CPs_standards_input <- CPs_standards |>
                    dplyr::select(Molecule, Batch_Name, RF) |>
                    tidyr::pivot_wider(names_from = Batch_Name, values_from = "RF") |>
                    dplyr::mutate(across(everything(), ~ tidyr::replace_na(., 0)))


                CPs_standards_sum_RF <- CPs_standards |>
                    dplyr::select(Batch_Name, Sum_RF_group) |>
                    dplyr::distinct() |>
                    dplyr::ungroup() |>
                    tidyr::pivot_wider(names_from = Batch_Name, values_from = "Sum_RF_group") |>
                    dplyr::mutate(across(everything(), ~ tidyr::replace_na(., 0)))


                # First populate combined_standard with CPs_standards_input
                combined_standard <- CPs_standards_input  |>
                    tibble::column_to_rownames(var = "Molecule") |>
                    as.matrix()



                # This performs deconvolution on all mixtures together. NEED TO CHECK IF perform_deconvolution ON SEPARATE SCCP, MCCP, LCCP is more correct!

                deconvolution <- combined_sample |>
                    #perform_deconvolution on only Relative_Area in the nested data frame
                    dplyr::mutate(result = purrr::map(data, ~ perform_deconvolution(dplyr::select(.x, Relative_Area), combined_standard, CPs_standards_sum_RF))) |>
                    dplyr::mutate(sum_Area = purrr::map_dbl(data, ~sum(.x$Area))) |>
                    dplyr::mutate(sum_deconv_RF = as.numeric(purrr::map(result, purrr::pluck("sum_deconv_RF")))) |>
                    dplyr::mutate(Concentration = sum_Area/sum_deconv_RF) |>
                    dplyr::mutate(deconv_coef = purrr::map(result, ~as_tibble(list(deconv_coef = .x$deconv_coef, Batch_Name = names(.x$deconv_coef))))) |>
                    dplyr::mutate(deconv_rsquared = as.numeric(purrr::map(result, purrr::pluck("deconv_rsquared")))) |>
                    dplyr::mutate(deconv_resolved = purrr::map(result, ~tibble::as_tibble(list(deconv_resolved = .x$deconv_resolved, Molecule = rownames(.x$deconv_resolved))))) |>
                    dplyr::select(-result)



                #### Calculate the concentration ####


                progress$set(value = 0.9, detail = "Calculating final results")



                # Store Samples_Concentration in the reactive value
                Samples_Concentration(deconvolution)

            }

            progress$set(value = 1, detail = "Complete")


            ### END: Deconvolution script



            output$downloadResults <- shiny::downloadHandler(
                filename = function() {
                    paste("CPquant_Results_", Sys.Date(), ".xlsx", sep = "")
                },
                content = function(file) {
                    # Create a new workbook
                    wb <- openxlsx::createWorkbook()

                    # Add worksheets
                    openxlsx::addWorksheet(wb, "Quantification")
                    openxlsx::addWorksheet(wb, "StandardsContribution")
                    openxlsx::addWorksheet(wb, "HomologueDistribution")

                    # Write data to worksheets
                    openxlsx::writeData(wb, "Quantification",
                                        deconvolution |>
                                            dplyr::select(Replicate_Name, Sample_Type, Concentration, deconv_rsquared) |>
                                            dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3)))
                    openxlsx::writeData(wb, "StandardsContribution",
                                        deconvolution |>
                                            unnest(deconv_coef) |>
                                            unnest_longer(c(deconv_coef, Batch_Name)) |>
                                            select(Replicate_Name, Batch_Name, deconv_coef) |>
                                            mutate(deconv_coef = deconv_coef * 100) |>
                                            pivot_wider(names_from = Batch_Name, values_from = deconv_coef))
                    openxlsx::writeData(wb, "HomologueDistribution",
                                        deconvolution |>
                                            mutate(data = map2(data, deconv_resolved, ~inner_join(.x, .y, by = "Molecule"))) |>
                                            mutate(data = map(data, ~ .x |> mutate(resolved_distribution = deconv_resolved / sum(deconv_resolved)))) |>
                                            select(-deconv_resolved) |>
                                            unnest(data) |>
                                            mutate(Deconvoluted_Distribution = as.numeric(resolved_distribution)) |>
                                            rename(Relative_Distribution = Relative_Area) |>
                                            select(-deconv_coef, -resolved_distribution, -Quantification_Group, -C_number, -Cl_number, -Area, -sum_Area,
                                                   -sum_deconv_RF, -Concentration, -deconv_resolved, -deconv_rsquared)
                    )


                    # Save workbook
                    openxlsx::saveWorkbook(wb, file)
                }
            )


            # Render table
            output$quantTable <- DT::renderDT({
                deconvolution |>
                    dplyr::select(Replicate_Name, Sample_Type, Concentration, deconv_rsquared) |> #remove to make compact df for pivot_wider
                    dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3)) |>
                    DT::datatable(
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

            # render sampleContributionPlot
            output$sampleContributionPlot <- renderPlotly({
                # How much contribution of each sample to the final deconvoluted homologue group pattern
                plot_data <- deconvolution |>
                    unnest(deconv_coef) |>
                    unnest_longer(c(deconv_coef, Batch_Name)) |>
                    select(Replicate_Name, Batch_Name, deconv_coef)

                # Create the plotly stacked bar plot
                plot_ly(plot_data,
                        x = ~Replicate_Name,
                        y = ~deconv_coef,
                        type = "bar",
                        color = ~Batch_Name,
                        colors = "Spectral") |>
                    layout(
                        title = list(
                            text = "Contributions from standards to deconvoluted pattern",
                            x = 0.5,  # Center the title
                            y = 0.95  # Position slightly down from top
                        ),
                        barmode = "stack",
                        xaxis = list(title = "Replicate Name"),
                        yaxis = list(title = "Relative Contribution",
                                     tickformat = ".2%"),
                        showlegend = TRUE,
                        legend = list(title = list(text = "Batch Name"))
                    )

            })




            ###########################################################QA/QC#######################################################

            QAQC <- deconvolution |>
                dplyr::select(Replicate_Name, Sample_Type, deconv_rsquared) |>
                dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3))

            if(input$calculateRecovery == "Yes") {
                # Recovery calculations
                recovery_data <- Skyline_output_filt |>
                    dplyr::filter(Isotope_Label_Type == "Quan",
                                  Molecule_List %in% c("IS", "RS"))


                RECOVERY <- recovery_data |>  # Calculate recovery
                    tidyr::pivot_wider(
                        id_cols = c(Replicate_Name, Sample_Type),
                        names_from = Molecule_List,
                        values_from = Area
                    ) |>
                    dplyr::mutate(across(c(IS, RS), ~replace_na(.x, 0)))

                # Calculate QC ratio
                qc_ratio <- RECOVERY |>
                    dplyr::filter(Sample_Type == "Quality Control") |>
                    dplyr::mutate(RatioStd = IS / RS) |>
                    dplyr::summarize(AverageRatio = mean(RatioStd, na.rm = TRUE))

                # Calculate sample recovery
                RECOVERY <- RECOVERY |>
                    dplyr::filter(Sample_Type %in% c("Unknown", "Blank")) |>
                    dplyr::mutate(
                        RatioSample = IS / RS,
                        Recovery = RatioSample / as.numeric(qc_ratio$AverageRatio),
                        RecoveryPercentage = round(Recovery * 100, 0)
                    ) |>
                    dplyr::select(Replicate_Name, Sample_Type, RecoveryPercentage)

                QAQC <- QAQC |>
                    dplyr::left_join(RECOVERY, by = c("Replicate_Name", "Sample_Type"))
            }

            # Render recovery table
            output$table_recovery <- DT::renderDT({
                DT::datatable(QAQC,
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



            if(input$calculateMDL == "Yes") {
                #MDL calculations (need to take into account if blank subtraction affect or not)
                if (input$blankSubtraction == "No"){

                    MDL_data <- deconvolution |>
                        dplyr::filter(Sample_Type == "Blank") |>
                        dplyr::summarize(
                            MDL_sumPCA = mean(Concentration) + 3 * sd(Concentration, na.rm = TRUE),
                            number_of_blanks = dplyr::n_distinct(Replicate_Name)
                        )
                }
                else if (input$blankSubtraction == "Yes, by avg area of blanks"){
                    DL_data <- deconvolution |>
                        dplyr::filter(Sample_Type == "Blank") |>
                        dplyr::summarize(
                            MDL_sumPCA = 3 * sd(Concentration, na.rm = TRUE),
                            number_of_blanks = dplyr::n_distinct(Replicate_Name)
                        )

                }
            }



            # Render MDL table
            output$table_MDL <- DT::renderDT({
                if (exists("MDL_data")) {
                    DT::datatable(MDL_data,
                                  options = list(
                                      pageLength = 10,
                                      dom = 't'
                                  ),
                                  rownames = FALSE)
                } else {
                    NULL
                }
            })








        })



        output$sampleSelectionUIOverlay <- renderUI({
            req(Samples_Concentration())

            # Get unique sample names
            sample_names <- unique(Samples_Concentration()$Replicate_Name)

            selectInput("selectedSamples", "Select samples to compare:",
                        choices = sample_names,
                        multiple = TRUE,
                        selected = NULL)
        })

        output$sampleSelectionUIComparisons <- renderUI({
            req(Samples_Concentration())

            # Get unique sample names
            sample_names <- unique(Samples_Concentration()$Replicate_Name)

            selectInput("selectedSamples", "Select samples to compare:",
                        choices = sample_names,
                        multiple = TRUE,
                        selected = NULL)
        })

        shiny::observeEvent(input$go2, {
            req(Samples_Concentration())  # Make sure the data exists

            withProgress(message = 'Generating plot...', value = 0, {

                Sample_distribution <- Samples_Concentration() |>
                    mutate(data = map2(data, deconv_resolved, ~inner_join(.x, .y, by = "Molecule"))) |>
                    mutate(data = map(data, ~ .x |> mutate(resolved_distribution = deconv_resolved / sum(deconv_resolved)))) |>
                    select(-deconv_resolved) |>
                    unnest(data)

                incProgress(0.5)

                if (input$plotHomologueGroups == "All Samples Overview") {
                    output$plotHomologuePatternStatic <- shiny::renderPlot({
                        #ggplot(Sample_distribution, aes(x = Molecule, y = resolved_distribution, fill = C_homologue)) +
                        ggplot(Sample_distribution, aes(x = Molecule, y = Relative_Area, fill = C_homologue)) +
                            geom_col() +
                            facet_wrap(~Replicate_Name) +
                            theme_minimal() +
                            theme(axis.text.x = element_blank()) +
                            labs(title = "Relative Distribution",
                                 x = "Homologue",
                                 y = "Relative Distribution")
                    })
                } else if (input$plotHomologueGroups == "Samples Overlay") {
                    output$plotHomologuePatternOverlay <- plotly::renderPlotly({
                        req(input$selectedSamples)
                        req(Sample_distribution)

                        # Filter for selected samples
                        selected_samples <- Sample_distribution |>
                            filter(Replicate_Name %in% input$selectedSamples) |>
                            mutate(Molecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)])))

                        req(nrow(selected_samples) > 0)

                        # Create a basic bar plot
                        p <- plot_ly(data = selected_samples,
                                     x = ~Molecule,
                                     #y = ~resolved_distribution,
                                     y = ~Relative_Area,
                                     color = ~Replicate_Name,
                                     type = "bar",
                                     text = ~paste(
                                         "Sample:", Replicate_Name,
                                         "<br>Homologue:", Molecule,
                                         "<br>Distribution:", round(Relative_Area, 4),
                                         "<br>C-atoms:", C_homologue
                                     ),
                                     hoverinfo = "text"
                        ) |>
                            layout(
                                title = "Sample Distribution Overlay",
                                xaxis = list(
                                    title = "Homologue",
                                    tickangle = 45
                                ),
                                yaxis = list(
                                    title = "Distribution"
                                ),
                                barmode = 'group',
                                showlegend = TRUE,
                                height = 600,
                                margin = list(b = 100)  # Add more bottom margin for rotated labels
                            )

                        p
                    })
                } else if (input$plotHomologueGroups == "Samples Panels") {
                    output$plotHomologuePatternComparisons <- plotly::renderPlotly({
                        req(input$selectedSamples)

                        # Filter data for selected samples and reshape data
                        selected_samples <- Sample_distribution |>
                            filter(Replicate_Name %in% input$selectedSamples) |>
                            mutate(Molecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)])))

                        # Get unique homologue groups for consistent coloring
                        homologue_groups <- unique(selected_samples$C_homologue)

                        # Create a list of plots, one for each Replicate_Name
                        plot_list <- selected_samples |>
                            split(selected_samples$Replicate_Name) |>
                            map(function(df) {
                                # Create base plot
                                p <- plot_ly()

                                # Add bars for each homologue group
                                for(hg in homologue_groups) {
                                    df_filtered <- df[df$C_homologue == hg,]
                                    p <- p |> add_trace(
                                        data = df_filtered,
                                        x = ~Molecule,
                                        #y = ~resolved_distribution,
                                        y = ~Relative_Area,
                                        name = hg,
                                        legendgroup = hg,
                                        showlegend = (df_filtered$Replicate_Name[1] == input$selectedSamples[1]),
                                        type = 'bar',
                                        opacity = 1
                                    )
                                }

                                # Add the black line for resolved_distribution
                                p <- p |> add_trace(
                                    data = df,
                                    x = ~Molecule,
                                    #y = ~Relative_Area,
                                    y = ~resolved_distribution,
                                    name = "Deconvoluted Distribution",
                                    #legendgroup = "RelativeArea",
                                    legendgroup = "DeconvDistr",
                                    showlegend = (df$Replicate_Name[1] == input$selectedSamples[1]),
                                    type = 'scatter',
                                    mode = 'lines+markers',
                                    line = list(color = 'black'),
                                    marker = list(color = 'black', size = 6),
                                    opacity = 0.7
                                )

                                # Add layout
                                p <- p |> layout(
                                    xaxis = list(
                                        title = "Homologue",
                                        tickangle = 45
                                    ),
                                    yaxis = list(title = "Value"),
                                    barmode = 'group',
                                    annotations = list(
                                        x = 0.5,
                                        y = 1.1,
                                        text = unique(df$Replicate_Name),
                                        xref = 'paper',
                                        yref = 'paper',
                                        showarrow = FALSE
                                    )
                                )

                                return(p)
                            })

                        # Combine the plots using subplot
                        subplot(plot_list,
                                nrows = ceiling(length(plot_list)/2),
                                shareX = TRUE,
                                shareY = TRUE) |>
                            layout(
                                title = "Sample Comparison",
                                showlegend = TRUE,
                                hovermode = 'closest',
                                hoverlabel = list(bgcolor = "white"),
                                barmode = 'group'
                            ) |>
                            config(displayModeBar = TRUE) |>
                            onRender("
                function(el) {
                    var plotDiv = document.getElementById(el.id);
                    plotDiv.on('plotly_legendclick', function(data) {
                        Plotly.restyle(plotDiv, {
                            visible: data.data[data.curveNumber].visible === 'legendonly' ? true : 'legendonly'
                        }, data.fullData.map((trace, i) => i).filter(i =>
                            data.fullData[i].legendgroup === data.fullData[data.curveNumber].legendgroup
                        ));
                        return false;
                    });
                }
            ")
                    })
                }



                incProgress(1)
            })
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
