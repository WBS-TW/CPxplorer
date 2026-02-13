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
#' @import stringr
#' @importFrom stats lm
#' @export


CPquant <- function(...){
    options(shiny.maxRequestSize = 500 * 1024^2)


    # =========================================================
    # UI
    # =========================================================


    ui <- shiny::navbarPage(
        "CPquant",
        shiny::tabPanel(
            "Quantification Inputs",
            shiny::fluidPage(
                shiny::sidebarLayout(
                    shiny::sidebarPanel(
                        width = 3,
                        shiny::fileInput("fileInput", "Import excel file from Skyline",
                                         accept = c('xlsx')),
                        shiny::textInput("quanUnit", "(Optional) Concentration unit:"),
                        shiny::radioButtons("quanSum",
                                            label = "Choose ions for quantification",
                                            choices = c("Quan only", "Sum Quan+Qual"),
                                            selected = "Quan only"),
                        shiny::radioButtons("blankSubtraction",
                                            label = "Subtraction with blank?",
                                            choices = c("Yes, by avg area of blanks", "No"),
                                            selected = "No"),
                        shiny::radioButtons("correctWithRS",
                                            label = "Correct with RS area?",
                                            choices = c("Yes", "No"),
                                            selected = "No"),
                        shiny::uiOutput("correctionUI"), # render UI if correctwithRS == "Yes"
                        shiny::radioButtons("calculateRecovery",
                                            label = "Calculate recovery? (req QC samples)",
                                            choices = c("Yes", "No"),
                                            selected = "No"),
                        shiny::radioButtons("calculateMDL",
                                            label = "Calculate MDL? (req blank samples)",
                                            choices = c("Yes", "No"),
                                            selected = "No"),
                        shiny::radioButtons("standardTypes",
                                            label = "Types of standards",
                                            choices = c("Group Mixtures"),
                                            selected = "Group Mixtures"),
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
            )
        ),
        shiny::tabPanel(
            "Input summary",
            shiny::sidebarPanel(
                width = 2,
                shiny::radioButtons("navSummary", "Choose tab:",
                                    choices = c("Std Calibration Curves",
                                                "Removed from Calibration",
                                                "Quan to Qual ratio",
                                                "Measured vs Theor Quan/Qual ratio"),
                                    selected = "Std Calibration Curves"),
                shiny::tags$hr(),
                shiny::tags$p("'Quan to Qual ratio' and 'Measured vs Theor Quan/Qual ratio' only works for Quan only option",
                              style = "font-style: italic;")
            ),
            shiny::mainPanel(
                width = 10,
                shiny::conditionalPanel(
                    condition = "input.navSummary == 'Std Calibration Curves'",
                    shiny::tags$h3("Standard calibration curves"),
                    plotly::plotlyOutput("CalibrationCurves")
                ),
                shiny::conditionalPanel(
                    condition = "input.navSummary == 'Removed from Calibration'",
                    shiny::tags$h3("Calibration series removed from quantification"),
                    DT::DTOutput("CalibrationRemoved")
                ),
                shiny::conditionalPanel(
                    condition = "input.navSummary == 'Quan to Qual ratio'",
                    shiny::tags$h3("Violin plots of Quant/Qual ions"),
                    plotly::plotlyOutput("RatioQuantToQual", height = "80vh", width = "100%")
                ),
                shiny::conditionalPanel(
                    condition = "input.navSummary == 'Measured vs Theor Quan/Qual ratio'",
                    shiny::tags$h3("Measured divided by Theoretical Quant/Qual ratios"),
                    plotly::plotlyOutput("MeasVSTheor", height = "80vh", width = "100%")
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
            shiny::fluidRow(
                shiny::column(
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
                shiny::column(
                    width = 10,
                    shiny::tags$div(
                        style = "margin-top: 20px;",
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
            )
        ),
        shiny::tabPanel(
            "QA/QC",
            shiny::mainPanel(
                DT::DTOutput("table_recovery"),
                shiny::br(),
                DT::DTOutput("table_MDL")
            )
        ),
        shiny::tabPanel(
            "Instructions",
            shiny::sidebarLayout(
                shiny::sidebarPanel(shiny::h3("Manual"), width = 3),
                shiny::mainPanel(
                    shiny::includeMarkdown(system.file("instructions_CPquant.md", package = "CPxplorer"))
                )
            )
        ),
        shiny::tabPanel(
            "Log",
            shiny::fluidPage(
                shiny::verbatimTextOutput("console_log"),
                shiny::tags$style("#console_log {max-height: 70vh;
                                    overflow-y: auto;
                                    background: #111827; /* dark */
                                    color: #e5e7eb;      /* light */
                                    font-family: monospace;
                                    font-size: 12px;
                                    padding: 12px;
                                    border-radius: 6px;}")
            )
        )
    )

    # =========================================================
    # Server
    # =========================================================

    server <- function(input, output, session) {

        # ---------------
        # In-memory console sink setup
        # ---------------
        log_buf <- character()
        tc <- textConnection("log_buf", open = "w", local = TRUE)
        sink(tc, append = TRUE)                     # stdout
        sink(tc, append = TRUE, type = "message")   # stderr (conditions)
        max_lines <- 600
        output$console_log <- shiny::renderText({
            shiny::invalidateLater(500, session)
            if (length(log_buf) > max_lines) {
                log_buf <<- tail(log_buf, max_lines)
            }
            paste(rev(log_buf), collapse = "\n")
        })
        # ---------------

        # ---- Quant unit passthrough ----
        quantUnit <- reactive({ input$quanUnit })

        # ---- 1) Load + minimal tidy (independent of removeSamples) ----
        base_df <- reactive({
            req(input$fileInput)

            df <- readxl::read_excel(
                input$fileInput$datapath,
                guess_max = 5000,
                na = c("", "NA", "#N/A", "N/A")
            )

            df <- df |>
                dplyr::rename(Replicate_Name = tidyr::any_of(c("Replicate Name", "ReplicateName"))) |>
                dplyr::rename(Sample_Type = tidyr::any_of(c("Sample Type", "SampleType"))) |>
                dplyr::rename(Molecule_List = tidyr::any_of(c("Molecule List", "MoleculeList"))) |>
                dplyr::rename(Mass_Error_PPM = tidyr::any_of(c("Mass Error PPM", "MassErrorPPM"))) |>
                dplyr::rename(Isotope_Label_Type = tidyr::any_of(c("Isotope Label Type", "IsotopeLabelType"))) |>
                dplyr::rename(Chromatogram_Precursor_MZ = tidyr::any_of(c("Chromatogram Precursor M/Z", "ChromatogramPrecursorMz"))) |>
                dplyr::rename(Analyte_Concentration = tidyr::any_of(c("Analyte Concentration", "AnalyteConcentration"))) |>
                dplyr::rename(Batch_Name = tidyr::any_of(c("Batch Name", "BatchName"))) |>
                dplyr::rename(Transition_Note = tidyr::any_of(c("Transition Note", "TransitionNote"))) |>
                dplyr::rename(Sample_Dilution_Factor = tidyr::any_of(c("Sample Dilution Factor", "SampleDilutionFactor"))) |>
                dplyr::mutate(
                    Analyte_Concentration = as.numeric(Analyte_Concentration),
                    Area = as.numeric(Area),
                    Area = tidyr::replace_na(Area, 0),
                    C_homologue = stringr::str_extract(Molecule, "C\\d+"),
                    Cl_homologue = stringr::str_extract(Molecule, "Cl\\d+"),
                    C_number = as.numeric(stringr::str_extract(C_homologue, "\\d+")),
                    Cl_number = as.numeric(stringr::str_extract(Cl_homologue, "\\d+")),
                    PCA = stringr::str_c(C_homologue, Cl_homologue, sep = ""),
                    Rel_Ab = as.numeric(purrr::map_chr(
                        Transition_Note, ~ {
                            matches <- stringr::str_match_all(.x, "\\{([^}]*)\\}")[[1]]
                            if (nrow(matches) >= 2) matches[2, 2] else NA_character_
                        }
                    ))
                )

            if (identical(input$standardTypes, "Group Mixtures")) {
                df <- df |>
                    dplyr::mutate(Quantification_Group = stringr::str_extract(Batch_Name, "^[^_]+"))
            }

            df
        })

        # ---- 2) Dynamic UI from base_df() ----
        output$defineVariables <- shiny::renderUI({
            req(base_df())
            defineVariablesUI(base_df())
        })

        output$correctionUI <- shiny::renderUI({
            req(base_df())
            if (input$correctWithRS == "Yes") {
                defineCorrectionUI(base_df())
            }
        })

        # ---- 3) Full processing with removal applied EARLY ----
        Skyline_output <- reactive({
            req(base_df())

            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message = "WAIT! Loading data...", value = 0)

            progress$set(value = 0.3, detail = "Starting from base data")
            df <- base_df()

            # Apply removal of samples BEFORE any downstream stats (RS, blanks, etc.)
            rmv <- input$removeSamples
            if (!is.null(rmv) && length(rmv) > 0) {
                df <- df |>
                    dplyr::filter(!Replicate_Name %in% rmv)
            }

            progress$set(value = 0.5, detail = "Applying quantification settings")

            # Sum all area of ions belong to same molecule if requested
            if (identical(input$quanSum, "Sum Quan+Qual")) {
                df <- df |>
                    dplyr::group_by(Replicate_Name, Molecule_List, Molecule) |>
                    dplyr::mutate(Area = sum(Area)) |>
                    dplyr::ungroup()
            }

            # RS normalization
            if (identical(input$correctWithRS, "Yes") && any(df$Molecule_List == "RS")) {
                if (!is.null(input$chooseRS) && nzchar(input$chooseRS)) {
                    df <- df |>
                        dplyr::group_by(Replicate_Name) |>
                        dplyr::mutate(
                            .rs_val = dplyr::first(Area[Molecule == input$chooseRS &
                                                            Molecule_List == "RS" &
                                                            Isotope_Label_Type == "Quan"]),
                            Area = dplyr::if_else(!is.na(.rs_val) & .rs_val > 0, Area / .rs_val, NA_real_)
                        ) |>
                        dplyr::ungroup() |>
                        dplyr::select(-.rs_val)
                }
            }

            # Blank subtraction — recompute after removal + RS correction
            if (identical(input$blankSubtraction, "Yes, by avg area of blanks")) {
                df_blank <- df |>
                    dplyr::filter(Sample_Type == "Blank") |>
                    dplyr::group_by(Molecule, Molecule_List, Isotope_Label_Type) |>
                    dplyr::summarise(AverageBlank = mean(Area, na.rm = TRUE), .groups = "drop") |>
                    dplyr::filter(!Molecule_List %in% c("IS", "RS", "VS"))

                df <- df |>
                    dplyr::left_join(df_blank,
                                     by = c("Molecule", "Molecule_List", "Isotope_Label_Type")) |>
                    dplyr::mutate(AverageBlank = tidyr::replace_na(AverageBlank, 0)) |>
                    dplyr::mutate(
                        Area = dplyr::case_when(
                            Sample_Type == "Unknown" ~ Area - AverageBlank,
                            TRUE ~ Area
                        ),
                        Area = ifelse(Area < 0, 0, Area)
                    )
            }

            if (identical(input$standardTypes, "Group Mixtures")) {
                df <- df |>
                    dplyr::mutate(
                        Quantification_Group = stringr::str_extract(Batch_Name, "^[^_]+")
                    )
            }

            progress$set(value = 1, detail = "Complete")
            print("Upload successful!")
            df
        })

        # ---- Event-gated input (for R^2) ----
        removeRsquared <- shiny::eventReactive(input$go, { as.numeric(input$removeRsquared) })
        Samples_Concentration <- reactiveVal() # store deconvolution for other tabs

        # ---- Raw table & plot (reflect auto removal) ----
        output$table_Skyline_output <- DT::renderDT({
            DT::datatable(Skyline_output(),
                          options = list(paging = TRUE, pageLength = 50))
        })

        output$plot_Skyline_output <- plotly::renderPlotly({
            plot_skyline_output(Skyline_output())
        })

        # =========================================================
        # Deconvolution pipeline — runs on Proceed button
        # =========================================================
        shiny::observeEvent(input$go, {

            progress <- shiny::Progress$new()
            on.exit(progress$close())

            progress$set(message = "WAIT! Processing data...", value = 0)

            # Skyline_output already reflects sample removals
            Skyline_output_filt <- Skyline_output()

            ##### PREPARE FOR DECONVOLUTION #######
            progress$set(value = 0.2, detail = "Preparing standards data")

            if (input$standardTypes == "Group Mixtures") {

                CPs_standards <- Skyline_output_filt |>
                    dplyr::filter(Sample_Type == "Standard",
                                  !Molecule_List %in% c("IS", "RS"),
                                  Isotope_Label_Type == "Quan",
                                  Batch_Name != "NA") |>
                    dplyr::mutate(
                        C_range = stringr::str_extract_all(Quantification_Group, "\\d+"),
                        C_min = as.numeric(purrr::map_chr(C_range, ~ .x[1])),
                        C_max = as.numeric(purrr::map_chr(C_range, function(x) { if (length(x) > 1) { x[2] } else { x[1] } }))
                    ) |>
                    dplyr::select(-C_range) |>
                    dplyr::group_by(Batch_Name, Sample_Type, Molecule, Molecule_List, C_number, Cl_number, PCA, Quantification_Group, C_min, C_max) |>
                    tidyr::nest() |>
                    dplyr::filter(C_number >= C_min & C_number <= C_max) |>
                    dplyr::mutate(models = purrr::map(data, ~ lm(Area ~ Analyte_Concentration, data = .x))) |>
                    dplyr::mutate(coef = purrr::map(models, coef)) |>
                    dplyr::mutate(RF = purrr::map_dbl(models, ~ coef(.x)["Analyte_Concentration"])) |>
                    dplyr::mutate(intercept = purrr::map(coef, purrr::pluck("(Intercept)"))) |>
                    dplyr::mutate(cal_rsquared = purrr::map(models, summary)) |>
                    dplyr::mutate(cal_rsquared = purrr::map(cal_rsquared, purrr::pluck("r.squared"))) |>
                    dplyr::select(-coef) |>
                    tidyr::unnest(c(RF, intercept, cal_rsquared)) |>
                    dplyr::mutate(RF = if_else(RF < 0, 0, RF)) |>
                    dplyr::mutate(cal_rsquared = ifelse(is.nan(cal_rsquared), 0, cal_rsquared)) |>
                    dplyr::mutate(RF = if_else(cal_rsquared < removeRsquared(), 0, RF)) |>
                    dplyr::ungroup() |>
                    dplyr::group_by(Batch_Name) |>
                    dplyr::mutate(Sum_RF_group = sum(RF, na.rm = TRUE)) |>
                    dplyr::ungroup()

                progress$set(value = 0.6, detail = "Preparing sample data")

                CPs_samples <- Skyline_output_filt |>
                    dplyr::filter(
                        Sample_Type %in% c("Unknown", "Blank"),
                        !Molecule_List %in% c("IS", "RS", "VS"),
                        Isotope_Label_Type == "Quan"
                    ) |>
                    dplyr::group_by(Replicate_Name) |>
                    dplyr::mutate(Relative_Area = Area / sum(Area, na.rm = TRUE)) |>
                    dplyr::ungroup() |>
                    dplyr::select(-Mass_Error_PPM, -Isotope_Label_Type, -Chromatogram_Precursor_MZ, -Analyte_Concentration, -Batch_Name) |>
                    dplyr::mutate(dplyr::across(Relative_Area, ~ replace(., is.nan(.), 0)))

                # Combined sample (nested) for per-sample deconvolution
                combined_sample <- CPs_samples |>
                    dplyr::group_by(Replicate_Name, Sample_Type) |>
                    tidyr::nest() |>
                    dplyr::ungroup()

                ###### Plots: Calibration Curves & Removed ######
                output$CalibrationCurves <- plotly::renderPlotly({
                    plot_calibration_curves(CPs_standards, quantUnit())
                })

                output$CalibrationRemoved <- DT::renderDT({
                    CPs_standards |>
                        dplyr::filter(RF <= 0) |>
                        dplyr::mutate(coef = purrr::map(models, coef)) |>
                        dplyr::mutate(RF = purrr::map_dbl(models, ~ coef(.x)["Analyte_Concentration"])) |>
                        dplyr::mutate(intercept = purrr::map(coef, purrr::pluck("(Intercept)"))) |>
                        dplyr::select(Batch_Name, Molecule, Quantification_Group, RF, intercept, cal_rsquared) |>
                        tidyr::unnest(c(RF, intercept)) |>
                        dplyr::mutate(across(where(is.numeric), ~ signif(.x, digits = 4))) |>
                        DT::datatable(options = list(pageLength = 40))
                })

                ###### Plots: Quan/Qual ratios ######
                output$RatioQuantToQual <- plotly::renderPlotly({
                    plot_quanqualratio(Skyline_output_filt)
                })
                output$MeasVSTheor <- plotly::renderPlotly({
                    plot_meas_vs_theor_ratio(Skyline_output_filt)
                })

                # ---- DECONVOLUTION ----
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

                # Matrix with molecules as rows, batches as columns
                combined_standard <- CPs_standards_input |>
                    tibble::column_to_rownames(var = "Molecule") |>
                    as.matrix()

                # This performs deconvolution on each sample
                deconvolution <- combined_sample |>
                    dplyr::mutate(result = purrr::map(
                        data,
                        ~ perform_deconvolution(dplyr::select(.x, Relative_Area), combined_standard, CPs_standards_sum_RF)
                    )) |>
                    dplyr::mutate(sum_Area = purrr::map_dbl(data, ~ sum(.x$Area))) |>
                    dplyr::mutate(sum_deconv_RF = as.numeric(purrr::map(result, purrr::pluck("sum_deconv_RF")))) |>
                    dplyr::mutate(Sample_Dilution_Factor = purrr::map_dbl(data, ~ dplyr::first(.x$Sample_Dilution_Factor))) |>
                    dplyr::mutate(
                        Concentration = dplyr::if_else(sum_deconv_RF > 0,
                                                       sum_Area / sum_deconv_RF * Sample_Dilution_Factor,
                                                       NA_real_)
                    ) |>
                    dplyr::mutate(Unit = quantUnit()) |>
                    dplyr::mutate(
                        deconv_coef = purrr::map(result, ~ tibble::tibble(
                            Batch_Name  = names(.x$deconv_coef),
                            deconv_coef = as.numeric(.x$deconv_coef)
                        )),
                        deconv_rsquared = as.numeric(purrr::map(result, purrr::pluck("deconv_rsquared"))),
                        deconv_resolved = purrr::map(result, ~ tibble::tibble(
                            Molecule         = names(.x$deconv_resolved),
                            deconv_resolved  = as.numeric(.x$deconv_resolved)
                        ))
                    ) |>
                    dplyr::select(-result)

                progress$set(value = 0.9, detail = "Calculating final results")

                # Store for downstream tabs
                Samples_Concentration(deconvolution)

                progress$set(value = 1, detail = "Complete")

                # ---- Download handler (Excel) ----
                output$downloadResults <- shiny::downloadHandler(
                    filename = function() {
                        paste("CPquant_Results_", Sys.Date(), ".xlsx", sep = "")
                    },
                    content = function(file) {
                        wb <- openxlsx::createWorkbook()

                        # Sheets
                        openxlsx::addWorksheet(wb, "Quantification")
                        openxlsx::addWorksheet(wb, "StandardsContribution")
                        openxlsx::addWorksheet(wb, "HomologueDistribution")

                        # Quantification
                        openxlsx::writeData(
                            wb, "Quantification",
                            deconvolution |>
                                dplyr::select(Replicate_Name, Sample_Type, Concentration, Unit, deconv_rsquared) |>
                                dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3))
                        )

                        # StandardsContribution
                        openxlsx::writeData(
                            wb, "StandardsContribution",
                            deconvolution |>
                                tidyr::unnest(deconv_coef) |>
                                dplyr::mutate(deconv_coef = deconv_coef * 100) |>
                                tidyr::pivot_wider(names_from = Batch_Name, values_from = deconv_coef)
                        )

                        # HomologueDistribution
                        openxlsx::writeData(
                            wb, "HomologueDistribution",
                            deconvolution |>
                                dplyr::mutate(data = purrr::map2(data, deconv_resolved, ~ dplyr::inner_join(.x, .y, by = "Molecule"))) |>
                                dplyr::mutate(data = purrr::map(data, ~ .x |>
                                                                    dplyr::mutate(resolved_distribution = deconv_resolved / sum(deconv_resolved)))) |>
                                dplyr::select(-deconv_resolved, -Sample_Dilution_Factor) |>
                                tidyr::unnest(data) |>
                                dplyr::mutate(Deconvoluted_Distribution = as.numeric(resolved_distribution)) |>
                                dplyr::rename(Relative_Distribution = Relative_Area) |>
                                dplyr::mutate(Molecule_Concentration = Deconvoluted_Distribution * Concentration) |>
                                dplyr::select(
                                    Replicate_Name, Sample_Type, Molecule_List, Molecule, C_homologue, Cl_homologue, PCA,
                                    Relative_Distribution, Deconvoluted_Distribution, Molecule_Concentration, Unit
                                )
                        )

                        openxlsx::saveWorkbook(wb, file)
                    }
                )

                # ---- Quant table ----
                output$quantTable <- DT::renderDT({
                    deconvolution |>
                        dplyr::select(Replicate_Name, Sample_Type, Concentration, Unit, deconv_rsquared) |>
                        dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3)) |>
                        DT::datatable(
                            filter = "top", extensions = c("Buttons", "Scroller"),
                            options = list(
                                scrollY = 650, scrollX = 500, deferRender = TRUE, scroller = TRUE,
                                buttons = list(
                                    list(extend = "excel", filename = "Samples_concentration", title = NULL,
                                         exportOptions = list(modifier = list(page = "all"))),
                                    list(extend = "csv", filename = "Samples_concentration", title = NULL,
                                         exportOptions = list(modifier = list(page = "all"))),
                                    list(extend = "colvis", targets = 0, visible = FALSE)
                                ),
                                dom = "lBfrtip",
                                fixedColumns = TRUE
                            ),
                            rownames = FALSE
                        )
                })

                # ---- Contribution plot ----
                output$sampleContributionPlot <- plotly::renderPlotly({
                    plot_sample_contribution(deconvolution)
                })

                # ---------------- QA/QC ----------------
                QAQC <- deconvolution |>
                    dplyr::select(Replicate_Name, Sample_Type, deconv_rsquared) |>
                    dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3))

                if (input$calculateRecovery == "Yes") {
                    recovery_data <- Skyline_output_filt |>
                        dplyr::filter(Isotope_Label_Type == "Quan",
                                      Molecule_List %in% c("IS", "RS"))

                    # Split IS and RS tables (keeping key context columns)
                    IS_tbl <- recovery_data %>%
                        dplyr::filter(Molecule_List == "IS") %>%
                        dplyr::select(Replicate_Name, Sample_Type, Molecule_IS = Molecule, Area_IS = Area)

                    RS_tbl <- recovery_data %>%
                        dplyr::filter(Molecule_List == "RS") %>%
                        dplyr::select(Replicate_Name, Sample_Type, Molecule_RS = Molecule, Area_RS = Area)

                    # Cross-join IS and RS within each replicate
                    ratio_all_pairs <- IS_tbl %>%
                        dplyr::inner_join(RS_tbl, by = c("Replicate_Name", "Sample_Type")) %>%
                        dplyr::mutate(
                            RatioISRS = dplyr::if_else(Area_RS > 0, Area_IS / Area_RS, NA_real_)
                        )

                    # Calculate QC ratio
                    qc_ratio <- ratio_all_pairs |>
                        dplyr::filter(Sample_Type == "Quality Control") |>
                        dplyr::group_by(Molecule_IS, Molecule_RS) |>
                        dplyr::summarize(AverageRatio = mean(RatioISRS, na.rm = TRUE), .groups = "drop")

                    ratio_with_sample <- ratio_all_pairs %>%
                        dplyr::filter(Sample_Type %in% c("Unknown", "Blank")) %>%
                        dplyr::mutate(
                            Molecule_IS = as.character(Molecule_IS),
                            Molecule_RS = as.character(Molecule_RS)
                        ) %>%
                        dplyr::left_join(
                            qc_ratio %>%
                                dplyr::mutate(
                                    Molecule_IS = as.character(Molecule_IS),
                                    Molecule_RS = as.character(Molecule_RS)
                                ),
                            by = c("Molecule_IS", "Molecule_RS")
                        ) %>%
                        dplyr::mutate(
                            Recovery = dplyr::if_else(!is.na(AverageRatio) & AverageRatio != 0,
                                                      round(RatioISRS / AverageRatio * 100, 1),
                                                      NA_real_)
                        ) %>%
                        dplyr::select(Replicate_Name, Sample_Type, Molecule_IS, Molecule_RS, Recovery)

                    QAQC <- QAQC |>
                        dplyr::left_join(ratio_with_sample, by = c("Replicate_Name", "Sample_Type"))
                }

                # Recovery table
                output$table_recovery <- DT::renderDT({
                    DT::datatable(
                        QAQC,
                        filter = "top",
                        extensions = c("Buttons", "Scroller"),
                        options = list(
                            scrollY = 650, scrollX = 500, deferRender = TRUE, scroller = TRUE,
                            buttons = list(
                                list(extend = "excel", filename = "Samples_recovery", title = NULL,
                                     exportOptions = list(modifier = list(page = "all"))),
                                list(extend = "csv", filename = "Samples_recovery", title = NULL,
                                     exportOptions = list(modifier = list(page = "all"))),
                                list(extend = "colvis", targets = 0, visible = FALSE)
                            ),
                            dom = "lBfrtip",
                            fixedColumns = TRUE
                        ),
                        rownames = FALSE
                    )
                })

                # MDL
                if (input$calculateMDL == "Yes") {
                    if (input$blankSubtraction == "No") {
                        MDL_data <- deconvolution |>
                            dplyr::filter(Sample_Type == "Blank") |>
                            dplyr::summarize(
                                MDL_sumPCA = mean(Concentration) + 3 * stats::sd(Concentration, na.rm = TRUE),
                                number_of_blanks = dplyr::n_distinct(Replicate_Name)
                            )
                    } else if (input$blankSubtraction == "Yes, by avg area of blanks") {
                        MDL_data <- deconvolution |>
                            dplyr::filter(Sample_Type == "Blank") |>
                            dplyr::summarize(
                                MDL_sumPCA = 3 * stats::sd(Concentration, na.rm = TRUE),
                                number_of_blanks = dplyr::n_distinct(Replicate_Name)
                            )
                    }
                }

                # MDL table
                output$table_MDL <- DT::renderDT({
                    if (exists("MDL_data")) {
                        DT::datatable(MDL_data, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
                    } else {
                        NULL
                    }
                })
            }
        })

        # ---- Homologue Group Patterns tab ----
        output$sampleSelectionUIOverlay <- renderUI({
            req(Samples_Concentration())
            sample_names <- unique(Samples_Concentration()$Replicate_Name)
            selectInput("selectedSamples", "Select samples to compare:",
                        choices = sample_names,
                        multiple = TRUE,
                        selected = NULL)
        })

        output$sampleSelectionUIComparisons <- renderUI({
            req(Samples_Concentration())
            sample_names <- unique(Samples_Concentration()$Replicate_Name)
            selectInput("selectedSamples", "Select samples to compare:",
                        choices = sample_names,
                        multiple = TRUE,
                        selected = NULL)
        })

        shiny::observeEvent(input$go2, {
            req(Samples_Concentration())

            withProgress(message = 'Generating plot...', value = 0, {
                Sample_distribution <- Samples_Concentration() |>
                    dplyr::mutate(data = purrr::map2(data, deconv_resolved, ~ dplyr::inner_join(.x, .y, by = "Molecule"))) |>
                    dplyr::mutate(data = purrr::map(data, ~ .x |>
                                                        dplyr::mutate(resolved_distribution = deconv_resolved / sum(deconv_resolved)))) |>
                    dplyr::select(-deconv_resolved, -Sample_Dilution_Factor) |>
                    tidyr::unnest(data) |>
                    dplyr::mutate(resolved_distribution = as.numeric(resolved_distribution)) |>
                    dplyr::mutate(deconv_resolved = as.numeric(deconv_resolved)) |>
                    dplyr::mutate(Molecule_concentration = resolved_distribution * Concentration)

                incProgress(0.5)

                if (input$plotHomologueGroups == "All Samples Overview") {
                    output$plotHomologuePatternStatic <- shiny::renderPlot({
                        ggplot2::ggplot(Sample_distribution, aes(x = Molecule, y = Relative_Area, fill = C_homologue)) +
                            ggplot2::geom_col() +
                            ggplot2::facet_wrap(~ Replicate_Name) +
                            ggplot2::theme_minimal() +
                            ggplot2::theme(axis.text.x = element_blank()) +
                            ggplot2::labs(
                                title = "Relative Distribution (before deconvolution)",
                                x = "Homologue",
                                y = "Relative Distribution"
                            )
                    })
                } else if (input$plotHomologueGroups == "Samples Overlay") {
                    output$plotHomologuePatternOverlay <- plotly::renderPlotly({
                        req(input$selectedSamples)
                        req(Sample_distribution)

                        selected_samples <- Sample_distribution |>
                            dplyr::filter(Replicate_Name %in% input$selectedSamples) |>
                            dplyr::mutate(Molecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)])))

                        req(nrow(selected_samples) > 0)

                        plotly::plot_ly(
                            data = selected_samples,
                            x = ~ Molecule,
                            y = ~ Relative_Area,   # switch to resolved_distribution if you prefer
                            color = ~ Replicate_Name,
                            type = "bar",
                            text = ~ paste(
                                "Sample:", Replicate_Name,
                                "<br>Homologue:", Molecule,
                                "<br>Distribution:", round(Relative_Area, 3),
                                "<br>C-atoms:", C_homologue
                            ),
                            hoverinfo = "text"
                        ) |>
                            layout(
                                xaxis = list(title = "Homologue", tickangle = 45),
                                yaxis = list(title = "Distribution"),
                                barmode = 'group',
                                showlegend = TRUE,
                                height = 600,
                                margin = list(b = 100)
                            )
                    })
                } else if (input$plotHomologueGroups == "Samples Panels") {
                    output$plotHomologuePatternComparisons <- plotly::renderPlotly({
                        req(input$selectedSamples)
                        plot_homologue_group_pattern_comparison(Sample_distribution, input$selectedSamples)
                    })
                }

                incProgress(1)
            })
        })

        # ---- Close the app when the session ends (non-interactive) ----
        if (!interactive()) {
            session$onSessionEnded(function() {
                stopApp()
                q("no")
            })
        }
    }

    # =========================================================
    # Run the application
    # =========================================================
    shinyApp(ui = ui, server = server)


}
