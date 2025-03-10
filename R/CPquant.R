#' CPquant: Shiny CPs Quantification for Skyline Output
#' @param ...
#'
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @import htmlwidgets
#' @rawNamespace import(ggplot2, except=c(last_plot))
#' @import readxl
#' @import nnls
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

    ui <- shiny::navbarPage("CPquant",
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
                                                    shiny::uiOutput("recoveryUI"), # render UI if correctwithRS == "Yes"
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
                                                        choices = c("Std Calibration Curves",
                                                                    "Removed from Calibration",
                                                                    "Quan to Qual ratio"),
                                                        selected = "Std Calibration Curves")
                                ),
                                shiny::mainPanel(
                                    width = 10,
                                    shiny::conditionalPanel(
                                        condition = "input.navSummary == 'Std Calibration Curves'",
                                        tags$h3("Standard calibration curves"),
                                        plotly::plotlyOutput("CalibrationCurves")
                                    ),
                                    shiny::conditionalPanel(
                                        condition = "input.navSummary == 'Removed from Calibration'",
                                        tags$h3("Calibration series removed from quantification"),
                                        DT::DTOutput("CalibrationRemoved")
                                    ),
                                    shiny::conditionalPanel(
                                        condition = "input.navSummary == 'Quan to Qual ratio'",
                                        tags$h3("Violin plots of Quant/Qual ions"),
                                        plotly::plotlyOutput("RatioQuantToQual", height = "80vh", width = "100%")
                                    )
                                )
                            ),

                            shiny::tabPanel(
                                "Quantification summary",
                                shiny::fluidPage(
                                    downloadButton("downloadResults", "Export all results to Excel"),
                                    shiny::tags$br(), shiny::tags$br(),
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

        # define RS
        output$recoveryUI <- shiny::renderUI({
            if(input$correctWithRS == "Yes") {
                defineRecoveryUI(Skyline_output()) }
        })
        #chooseRS <- shiny::reactive(as.character(input$chooseRS))

        Skyline_output <- reactive({
            req(input$fileInput) #requires that the input is available

            # Create a Progress object
            progress <- shiny::Progress$new()
            on.exit(progress$close())

            progress$set(message = "WAIT! Loading data...", value = 0)

            # Read the Skyline output Excel file
            progress$set(value = 0.3, detail = "Reading Excel file")
            df <- readxl::read_excel(input$fileInput$datapath) #outputs a tibble

            progress$set(value = 0.6, detail = "Processing data")
            df <- df |>
                dplyr::rename(Replicate_Name = tidyr::any_of(c("Replicate Name", "ReplicateName"))) |>
                dplyr::rename(Sample_Type = tidyr::any_of(c("Sample Type", "SampleType"))) |>
                dplyr::rename(Molecule_List = tidyr::any_of(c("Molecule List", "MoleculeList"))) |>
                dplyr::rename(Mass_Error_PPM = tidyr::any_of(c("Mass Error PPM", "MassErrorPPM"))) |>
                dplyr::rename(Isotope_Label_Type = tidyr::any_of(c("Isotope Label Type", "IsotopeLabelType"))) |>
                dplyr::rename(Chromatogram_Precursor_MZ = tidyr::any_of(c("Chromatogram Precursor M/Z", "ChromatogramPrecursorMz"))) |>
                dplyr::rename(Analyte_Concentration = tidyr::any_of(c("Analyte Concentration", "AnalyteConcentration"))) |>
                dplyr::rename(Batch_Name = tidyr::any_of(c("Batch Name", "BatchName"))) |>
                dplyr::mutate(Analyte_Concentration = as.numeric(Analyte_Concentration)) |>
                dplyr::mutate(Area = as.numeric(Area)) |>
                dplyr::mutate(Area = replace_na(Area, 0)) |>
                dplyr::mutate(C_homologue = stringr::str_extract(Molecule, "C\\d+"),
                              Cl_homologue = stringr::str_extract(Molecule, "Cl\\d+"),
                              C_number = as.numeric(stringr::str_extract(C_homologue, "\\d+")),
                              Cl_number = as.numeric(stringr::str_extract(Cl_homologue, "\\d+")),
                              PCA = stringr::str_c(C_homologue, Cl_homologue, sep = ""))

            progress$set(value = 0.8, detail = "Applying corrections")


            # Normalize data based on 'Correct with RS' input
            if (input$correctWithRS == "Yes" && any(df$Molecule_List == "RS")) {
                # Only proceed with RS correction if chooseRS input is available
                if (!is.null(input$chooseRS) && input$chooseRS != "") {
                    # Use input$chooseRS directly instead of the reactive
                    df <- df |>
                        dplyr::group_by(Replicate_Name) |>
                        dplyr::mutate(Area = Area / first(Area[Molecule == input$chooseRS &
                                                                   Molecule_List == "RS" &
                                                                   Isotope_Label_Type == "Quan"])) |>
                        dplyr::ungroup()
                }
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
                    dplyr::mutate(Quantification_Group = stringr::str_extract(Batch_Name, "^[^_]+")) #extract the initial sequence of characters up to, but not including, the first underscore.
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
        removeSamples <- shiny::eventReactive(input$go, {as.character(input$removeSamples)})
        #chooseRS <- shiny::eventReactive(input$go, {as.character(input$chooseRS)})
        Samples_Concentration <- reactiveVal() # Create a reactive value to store deconvolution object into Samples_Concentration() to allow other to access after observeEvent.


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
        # plots.R function
        output$plot_Skyline_output <- plotly::renderPlotly({
            plot_skyline_output(Skyline_output())

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
                                  Isotope_Label_Type == "Quan", # use only Quan ions
                                  Batch_Name != "NA") |>
                    dplyr::mutate(
                        C_range = stringr::str_extract_all(Quantification_Group, "\\d+"), #extract all sequences of digits
                        C_min = as.numeric(purrr::map_chr(C_range, ~.x[1])),
                        C_max = as.numeric(purrr::map_chr(C_range, function(x) {if(length(x) > 1) {x[2]} else {x[1]}}))) |>
                    dplyr::select(-C_range) |>
                    dplyr::group_by(Batch_Name, Sample_Type, Molecule, Molecule_List, C_number, Cl_number, PCA, Quantification_Group, C_min, C_max) |>
                    tidyr::nest() |>
                    dplyr::filter(C_number >= C_min & C_number <= C_max) |> # make sure C_number stays within the Quantification_Group chain length
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
                    dplyr::mutate(RF = if_else(rsquared < removeRsquared(), 0, RF)) |> #replace RF with 0 if rsquared is below removeRsquared()
                    dplyr::ungroup() |>
                    dplyr::group_by(Batch_Name) |> #grouping by the standards
                    dplyr::mutate(Sum_RF_group = sum(RF, na.rm = TRUE)) |>
                    dplyr::ungroup()



                # Prepare for deconvolution of samples
                progress$set(value = 0.6, detail = "Preparing sample data")



                CPs_samples <- Skyline_output_filt |>
                    dplyr::filter(
                        #Sample_Type == "Unknown",
                        Sample_Type %in% c("Unknown", "Blank"), #include both unknown and blank
                        !Molecule_List %in% c("IS", "RS", "VS"), # remove IS, RS, VS
                        Isotope_Label_Type == "Quan") |>
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


                ###### Plot calibration curves ######
                # plots.R function
                output$CalibrationCurves <- plotly::renderPlotly({
                    plot_calibration_curves(CPs_standards)
                })

                ###### Plot CalibrationRemoved ######
                output$CalibrationRemoved <- DT::renderDT({
                    CPs_standards |>
                        dplyr::filter(RF <= 0) |>
                        dplyr::mutate(coef = purrr::map(models, coef)) |>
                        dplyr::mutate(RF = purrr::map_dbl(models, ~ coef(.x)["Analyte_Concentration"]))|> #get the slope which will be the RF
                        dplyr::mutate(intercept = purrr::map(coef, purrr::pluck("(Intercept)"))) |>
                        dplyr::select(Batch_Name, Molecule, Quantification_Group,RF, intercept,  rsquared) |>
                        tidyr::unnest(c(RF, intercept)) |>
                        mutate(across(where(is.numeric), ~ signif(.x, digits = 4))) |>
                        DT::datatable(options = list(pageLength = 40))
                })


                ###### Plot Quan/Qual ratios ######
                #plots.R function
                output$RatioQuantToQual <- plotly::renderPlotly({
                    plot_quanqualratio(Skyline_output_filt)
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



                # This performs deconvolution on all mixtures together.

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


            # download results from deconvolution to different excel sheets
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
                    dplyr::select(Replicate_Name, Sample_Type, Concentration, deconv_rsquared) |> #select to make compact df for pivot_wider
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
               plot_sample_contribution(deconvolution)})



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
                        ggplot(Sample_distribution, aes(x = Molecule, y = Relative_Area, fill = C_homologue)) +
                            geom_col() +
                            facet_wrap(~Replicate_Name) +
                            theme_minimal() +
                            theme(axis.text.x = element_blank()) +
                            labs(title = "Relative Distribution (before deconvolution)",
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
                                         "<br>Distribution:", round(Relative_Area, 3),
                                         "<br>C-atoms:", C_homologue
                                     ),
                                     hoverinfo = "text"
                        ) |>
                            layout(
                                #title = "Sample Distribution Overlay",
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
                        #plots.R function
                        plot_homologue_group_pattern_comparison(Sample_distribution, input$selectedSamples)
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
