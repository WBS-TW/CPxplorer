#' CPions
#'
#' @param ...
#'
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @import shinythemes
#' @import DT
#' @import plotly
#' @import dplyr
#' @rawNamespace import(ggplot2, except=c(last_plot))
#' @import purrr
#' @import readr
#' @import stringr
#' @import tibble
#' @import tidyr
#' @import readxl
#' @import enviPat
#' @import markdown
#' @export


CPions <- function(...){

table_download_controls <- function(csv_id, xlsx_id) {
    shiny::div(
        style = "margin-bottom: 10px; display: flex; gap: 8px; flex-wrap: wrap;",
        shiny::downloadButton(csv_id, "Download CSV"),
        shiny::downloadButton(xlsx_id, "Download Excel")
    )
}

#--------------------------------UI function----------------------------------#
ui <- shiny::navbarPage(
    "CPions",
    theme = shinythemes::shinytheme('spacelab'),
    shiny::tabPanel("Normal settings",
                    shiny::fluidPage(shiny::sidebarLayout(
                        shiny::sidebarPanel(
                            shiny::numericInput("Cmin", "C atoms min (allowed 3-40)", value = 10, min = 3, max = 40),
                            shiny::numericInput("Cmax", "C atoms max (allowed 4-40)", value = 30, min = 4, max = 40),
                            shiny::numericInput("Clmin", "Cl atoms min (allowed 1-15))", value = 3, min = 1, max = 15),
                            shiny::numericInput("Clmax", "Cl atoms max (allowed 1-15)", value = 15, min = 1, max = 15),
                            shiny::numericInput("Brmin", "Br atoms min (allowed 1-15))", value = 1, min = 1, max = 15),
                            shiny::numericInput("Brmax", "Br atoms max (allowed 1-15)", value = 4, min = 1, max = 15),
                            shiny::br(),
                            shiny::selectInput("Adducts", "Add adducts/fragments",
                                        choices = c("[PCA-Cl]-",
                                                    "[PCA-H]-",
                                                    "[PCA-HCl]-",
                                                    "[PCA-Cl-HCl]-",
                                                    "[PCA-2Cl-HCl]-",
                                                    "[PCA+Cl]-",
                                                    "[PCO-Cl]-",
                                                    "[PCO-HCl]-",
                                                    "[PCO-H]-",
                                                    "[PCO+Cl]-",
                                                    "[PCA+Br]-",
                                                    "[BCA+Cl]-",
                                                    "[BCA-Cl]-",
                                                    "[PCA-Cl-HCl]+",
                                                    "[PCA-Cl-2HCl]+",
                                                    "[PCA-Cl-3HCl]+",
                                                    "[PCA-Cl-4HCl]+"
                                        ),
                                        selected = "[PCA+Cl]-",
                                        multiple = TRUE,
                                        selectize = TRUE,
                                        width = NULL,
                                        size = NULL),
                            shiny::numericInput("threshold", "Isotope rel ab threshold (1-99%)", value = 5, min = 1, max = 99),
                            shiny::textAreaInput("ISRS_input", "Optional: add ion formula for IS/RS",
                                               placeholder = "Input the [M+adduct] ion formula. See Instructions", height = "150px"),
                            shiny::actionButton("go1", "Submit", width = "100%"),
                            width = 3),
                        shiny::mainPanel(
                            table_download_controls("download_norm_csv", "download_norm_xlsx"),
                            DT::DTOutput("Table_norm", width = "100%")
                        )
                    )
                    )),
    shiny::tabPanel("Advanced settings",
                    shiny::fluidPage(shiny::sidebarLayout(
                        shiny::sidebarPanel(
                            shiny::numericInput("Cmin_adv", "C atoms min (allowed 3-40)", value = 10, min = 3, max = 40),
                            shiny::numericInput("Cmax_adv", "C atoms max (allowed 4-40)", value = 30, min = 4, max = 40),
                            shiny::numericInput("Clmin_adv", "Cl atoms min (allowed 1-15))", value = 3, min = 1, max = 15),
                            shiny::numericInput("Clmax_adv", "Cl atoms max (allowed 1-15)", value = 15, min = 1, max = 15),
                            shiny::numericInput("Brmin_adv", "Br atoms min (allowed 1-15))", value = 1, min = 1, max = 15),
                            shiny::numericInput("Brmax_adv", "Br atoms max (allowed 1-15)", value = 4, min = 1, max = 15),
                            shiny::br(),
                            shiny::selectInput("Compclass_adv", "Compound Class",
                                        choices = c("PCA", "PCO","BCA"),
                                        selected = "PCA",
                                        multiple = TRUE,
                                        selectize = TRUE,
                                        width = NULL,
                                        size = NULL),
                            shiny::selectInput("Adducts_adv", "Which adduct",
                                        choices = c("-Cl", "-H", "-HCl", "-Cl-HCl","-Cl-2HCl", "-Cl-3HCl", "-2Cl-HCl", "+Cl","+Br", "+F"),
                                        selected = "+Cl",
                                        multiple = TRUE,
                                        selectize = TRUE,
                                        width = NULL,
                                        size = NULL),
                            shiny::selectInput("Charge_adv", "Which charge",
                                        choices = c("-", "+"),
                                        selected = "-",
                                        multiple = FALSE,
                                        selectize = TRUE,
                                        width = NULL,
                                        size = NULL),
                            shiny::selectInput("TP_adv", "Transformation product",
                                        choices = c("None", "-Cl+OH", "-2Cl+2OH", "-H+OH", "-2H+2OH", "-2H+O", "-H+SO4H"),
                                        selected = "None",
                                        multiple = TRUE,
                                        selectize = TRUE,
                                        width = NULL,
                                        size = NULL),
                            shiny::numericInput("threshold_adv", "Isotope rel ab threshold (0-99%)", value = 5, min = 0, max = 99),
                            shiny::textAreaInput("ISRS_input_adv", "Optional: add ion formula for IS/RS",
                                                 placeholder = "Input the [M+adduct] ion formula. See Instructions" , height = "150px"),

                            shiny::actionButton("go_adv", "Submit", width = "100%"),
                            width = 3),
                        shiny::mainPanel(
                            table_download_controls("download_adv_csv", "download_adv_xlsx"),
                            DT::DTOutput("Table_adv", width = "100%")
                        )
                    )
                    )),
    shiny::tabPanel("Interfering ions",
                    shiny::fluidPage(shiny::sidebarLayout(
                        shiny::sidebarPanel(
                            shiny::numericInput("MSresolution", "MS Resolution", value = 20000, min = 100, max = 5000000),
                            shiny::radioButtons("interfere_mode", label = "From Normal or Advanced settings", choices = c("normal", "advanced"), selected = "normal"),
                            shiny::actionButton("go2", "Calculate", width = "100%"),
                            width = 2
                        ),
                        shiny::mainPanel(
                            plotly::plotlyOutput("Plotly"),
                            plotly::plotlyOutput("Plotly2"),
                            table_download_controls("download_interference_csv", "download_interference_xlsx"),
                            DT::DTOutput("Table2", width = "100%")

                        )
                    )
                    )),
    shiny::tabPanel("Skyline",
                    shiny::fluidPage(shiny::sidebarLayout(
                        shiny::sidebarPanel(
                            #shiny::radioButtons("skylineoutput", label = "Output table", choices = c("mz", "IonFormula")),
                            shiny::radioButtons("skylineoutput", label = "Output table", choices = c("mz")),
                            shiny::radioButtons("skyline_mode", label = "From Normal or Advanced settings", choices = c("normal", "advanced"), selected = "normal"),
                             shiny::radioButtons("QuantIon", label = "Use as Quant Ion", choices = c("Most intense", "Interference-filtered"),
                                                 selected = "Most intense"),
                             shiny::conditionalPanel(
                                 condition = "input.QuantIon == 'Interference-filtered'",
                                 shiny::helpText("When using interference filtering, choose a strategy to select the Quant ion and set the preferred number of Qual ions."),
                                 shiny::selectInput(
                                     "skyline_strategy",
                                     label = "Selection strategy",
                                     choices = c(
                                        "Highest abundance Quan" = "abundance",
                                        "Balanced" = "balanced",
                                        "Least interference Quan" = "interference"
                                     ),
                                     selected = "abundance"
                                 ),
                                 shiny::numericInput(
                                     "skyline_qual_n",
                                     label = "Preferred number of Qual ions",
                                     value = 2,
                                     min = 0,
                                     max = 20,
                                     step = 1
                                 )
                             ),
                            shiny::actionButton("go3", "Transition List", width = "100%"),
                            width = 4
                        ),
                        shiny::mainPanel(
                            shiny::div(
                                style = "margin-bottom: 10px; display: flex; gap: 8px; flex-wrap: wrap;",
                                shiny::downloadButton("download_iontable_csv", "Download Full Ion Table CSV"),
                                shiny::downloadButton("download_iontable_xlsx", "Download Full Ion Table Excel")
                            ),
                            shiny::div(
                                style = "margin-bottom: 10px; display: flex; gap: 8px; flex-wrap: wrap;",
                                shiny::downloadButton("download_skyline_csv", "Download Skyline Transition List CSV"),
                                shiny::downloadButton("download_skyline_xlsx", "Download Skyline Transition List Excel")
                            ),
                            DT::DTOutput("Table3", width = "100%")

                        )
                    )
                    )),
    shiny::tabPanel(
        "Instructions",
        shiny::sidebarLayout(
            shiny::sidebarPanel(shiny::h3("Manual"),
                                width = 3),
            shiny::mainPanel(
                shiny::includeMarkdown(system.file("instructions_CPions.md", package = "CPxplorer"))
            )
        )
    )
)



#---------------------------Shiny Server function------------------------------#



server = function(input, output, session) {
    CPions_server(input, output, session)
    # Close the app when the session ends
    if(!interactive()) {
        session$onSessionEnded(function() {
            stopApp()
            q("no")
        })
    }

}

shiny::shinyApp(ui, server)

}

CPions_server <- function(input, output, session) {
    suppress_dt_size_warning <- function(expr) {
        old_options <- options(DT.warn.size = FALSE)
        on.exit(options(old_options), add = TRUE)
        force(expr)
    }

    register_table_downloads <- function(data_reactive, csv_output_id, xlsx_output_id, file_stub) {
        output[[csv_output_id]] <- shiny::downloadHandler(
            filename = function() {
                paste0(file_stub, ".csv")
            },
            content = function(file) {
                data <- data_reactive()
                shiny::req(data)
                readr::write_csv(data, file)
            }
        )

        output[[xlsx_output_id]] <- shiny::downloadHandler(
            filename = function() {
                paste0(file_stub, ".xlsx")
            },
            content = function(file) {
                data <- data_reactive()
                shiny::req(data)
                openxlsx::write.xlsx(data, file = file, overwrite = TRUE)
            }
        )
    }

    # Set reactive values from user input

    # GENERAL
    MSresolution <- shiny::eventReactive(input$go2, {as.integer(input$MSresolution)})
    CP_allions_compl2 <- shiny::reactiveVal(NULL) # save as global object after calculate the interfering ions for skyline
    CP_allions_skyline_reactive <- shiny::reactiveVal(NULL)

    # NORMAL
    selectedAdducts <- shiny::eventReactive(input$go1, {as.character(input$Adducts)})

    C <- shiny::eventReactive(input$go1, {as.integer(input$Cmin:input$Cmax)})
    Cl <- shiny::eventReactive(input$go1, {as.integer(input$Clmin:input$Clmax)})
    Clmin <- shiny::eventReactive(input$go1, {as.integer(input$Clmin)})
    Clmax <- shiny::eventReactive(input$go1, {as.integer(input$Clmax)})
    Br <- shiny::eventReactive(input$go1, {as.integer(input$Brmin:input$Brmax)})
    Brmin <- shiny::eventReactive(input$go1, {as.integer(input$Brmin)})
    Brmax <- shiny::eventReactive(input$go1, {as.integer(input$Brmax)})
    threshold <- shiny::eventReactive(input$go1, {as.integer(input$threshold)})
    ISRS_input <- shiny::eventReactive(input$go1, {as.character(input$ISRS_input)})


    # ADVANCED
    selectedClass_adv <- shiny::eventReactive(input$go_adv, {as.character((input$Compclass_adv))})
    selectedAdducts_adv <- shiny::eventReactive(input$go_adv, {as.character(input$Adducts_adv)})
    selectedCharge_adv <- shiny::eventReactive(input$go_adv, {as.character((input$Charge_adv))})
    selectedTP_adv <- shiny::eventReactive(input$go_adv, {as.character((input$TP_adv))})

    C_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Cmin_adv:input$Cmax_adv)})
    Cl_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Clmin_adv:input$Clmax_adv)})
    Clmin_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Clmin_adv)})
    Clmax_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Clmax_adv)})
    Br_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Brmin_adv:input$Brmax_adv)})
    Brmin_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Brmin_adv)})
    Brmax_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Brmax_adv)})
    threshold_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$threshold_adv)})
    ISRS_input_adv <- shiny::eventReactive(input$go_adv, {as.character(input$ISRS_input_adv)})




    #----Outputs_Start
    # NORMAL
    CP_allions_glob <- shiny::eventReactive(input$go1, {

        # Create a Progress bar object
        progress <- shiny::Progress$new()

        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        progress$set(message = "Calculating", value = 0)

        Adducts <- as.character(selectedAdducts())

        # function to get adducts or fragments
        CP_allions_template <- data.frame(Molecule_Formula = character(), Halo_perc = double())
        CP_allions_inputs <- vector("list", length(Adducts))
        for (i in seq_along(Adducts)) {
            progress$inc(1/length(Adducts), detail = paste0("Adduct: ", Adducts[i], " . Please wait.."))
            if(stringr::str_detect(Adducts[i], "\\bBCA\\b")){
                CP_allions_inputs[[i]] <- getAdduct_BCA(adduct_ions = Adducts[i], C = C(), Cl = Cl(), Clmax = Clmax(),
                                       Br = Br(), Brmax = Brmax(), threshold = threshold())
            } else {
                CP_allions_inputs[[i]] <- getAdduct_normal(adduct_ions = Adducts[i], C = C(), Cl = Cl(), Clmax = Clmax(), threshold = threshold())
            }
        }

        CP_allions <- combine_cpions_tables(CP_allions_inputs, CP_allions_template)


        # Add ISRS if textinput is not empty ""
        if(ISRS_input() != ""){
        CP_allions <- addISRS(ISRS_input(), CP_allions, threshold())
        }


        return(CP_allions)
    })

    shiny::observeEvent(input$go1, {
        output$Table_norm <- DT::renderDT({
            suppress_dt_size_warning(
                DT::datatable(CP_allions_glob(),
                              filter = "top", extensions = c("Buttons", "Scroller"),
                              options = list(scrollY = 650,
                                             scrollX = 500,
                                             deferRender = TRUE,
                                             scroller = TRUE,
                                             buttons = list(list(extend = "colvis", targets = 0, visible = FALSE)),
                                             dom = "lBfrtip",
                                             fixedColumns = TRUE),
                              rownames = FALSE)
            )
        }, server = TRUE)
    })
    # go1 end
    register_table_downloads(CP_allions_glob, "download_norm_csv", "download_norm_xlsx", "CPions_normal_settings")


    # ADVANCED
    CP_allions_glob_adv <- shiny::eventReactive(input$go_adv, {

        # Create a Progress bar object
        progress <- shiny::Progress$new()

        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        progress$set(message = "Calculating", value = 0)

        Class <- as.character(selectedClass_adv())
        Adducts <- as.character(selectedAdducts_adv())
        Charge <- as.character(selectedCharge_adv())
        TP <- as.character(selectedTP_adv())

        # function to get adducts or fragments
        CP_allions_template <- data.frame(Molecule_Formula = character(), Molecule_Halo_perc = double())
        CP_allions_inputs <- vector("list", length(Class) * length(Adducts) * length(TP))
        input_idx <- 0L

        # nested for loop to get all combinations of Class, Adducts, TP
        for (i in seq_along(Class)) {
            progress$inc(1/length(Class), detail = paste0("Compound Class: ", Class[i], " . Please wait.."))

            for (j in seq_along(Adducts)) {

                for (k in seq_along(TP)) {
                    input_idx <- input_idx + 1L
                    CP_allions_inputs[[input_idx]] <- getAdduct_advanced(Class = Class[i], Adduct_Ion = Adducts[j], TP = TP[k], Charge = Charge,
                                                C = C_adv(), Cl = Cl_adv(), Clmax = Clmax_adv(), Br = Br_adv(), Brmax = Brmax_adv(),
                                                threshold = threshold_adv())
                }
            }
        }

        CP_allions <- combine_cpions_tables(CP_allions_inputs, CP_allions_template)

        # Add ISRS if textinput is not empty ""
        if(ISRS_input_adv() != ""){
            CP_allions <- addISRS(ISRS_input_adv(), CP_allions, threshold_adv())
        }

        return(CP_allions)

    })

    ### go_adv: Calculate the isotopes from initial settings tab ###
    shiny::observeEvent(input$go_adv, {
        output$Table_adv <- DT::renderDT({
            suppress_dt_size_warning(
                DT::datatable(CP_allions_glob_adv(),
                              filter = "top", extensions = c("Buttons", "Scroller"),
                              options = list(scrollY = 650,
                                             scrollX = 500,
                                             deferRender = TRUE,
                                             scroller = TRUE,
                                             buttons = list(list(extend = "colvis", targets = 0, visible = FALSE)),
                                             dom = "lBfrtip",
                                             fixedColumns = TRUE),
                              rownames = FALSE)
            )
        }, server = TRUE)
    })
    # go_adv end
    register_table_downloads(CP_allions_glob_adv, "download_adv_csv", "download_adv_xlsx", "CPions_advanced_settings")

    ##################################################################
    ############ go2: Calculates the interfering ions tab ############
    ##################################################################
    shiny::observeEvent(input$go2, {


        CP_allions_interfere <- if (input$interfere_mode == "normal") {
            CP_allions_glob()
        } else {
            CP_allions_glob_adv()
        }


        CP_allions_interfere <- compute_interference(CP_allions_interfere, MSresolution())

        # populates CP_allions_compl2 so it can be use for skyline tab
        CP_allions_compl2(CP_allions_interfere)

        # Output scatterplot: #Cl vs #C  if Br exists
        if ("79Br" %in% names(CP_allions_compl2) == TRUE){
            output$Plotly <- plotly::renderPlotly(
                p <- CP_allions_interfere |>
                    dplyr::mutate(`79Br` = tidyr::replace_na(`79Br`, 0)) |>
                    dplyr::mutate(`81Br` = tidyr::replace_na(`81Br`, 0)) |>
                    plotly::plot_ly(
                        x = ~ (`12C`+`13C`),
                        y = ~(`35Cl`+`37Cl`+`79Br`+`81Br`),
                        type = "scatter",
                        mode = "markers",
                        color = ~interference,
                        hoverinfo = "text",
                        hovertext = paste("Molecule_Formula:", CP_allions_interfere$Molecule_Formula,
                                          '<br>',
                                          "Adduct/Fragment ion:", paste0(CP_allions_interfere$Adduct, CP_allions_interfere$Isotopologue),
                                          '<br>',
                                          "Ion Formula:", CP_allions_interfere$Adduct_Formula,
                                          '<br>',
                                          "Adduct isotopes:", paste0("[12C]:", CP_allions_interfere$`12C`, "  [13C]:", CP_allions_interfere$`13C`,
                                                                     "  [35Cl]:", CP_allions_interfere$`35Cl`, "  [37Cl]:", CP_allions_interfere$`37Cl`, " [79Br]:", CP_allions_interfere$`79Br`, " [81Br]:", CP_allions_interfere$`81Br`))
                    )
                |>
                    plotly::layout(xaxis = list(title = "Number of carbons (12C+13C)"),
                                   yaxis = list(title = "Number of halogens (35Cl+37Cl+79Br+81Br)"),
                                   legend=list(title=list(text='<b> Interference at MS res? </b>')))
            )
        } else { #if there are no bromines
            output$Plotly <- plotly::renderPlotly(
                p <- CP_allions_interfere |>
                    plotly::plot_ly(
                        x = ~ (`12C`+`13C`),
                        y = ~(`35Cl`+`37Cl`),
                        type = "scatter",
                        mode = "markers",
                        color = ~interference,
                        hoverinfo = "text",
                        hovertext = paste("Molecule_Formula:", CP_allions_interfere$Molecule_Formula,
                                          '<br>',
                                          "Adduct/Fragment ion:", paste0(CP_allions_interfere$Adduct, CP_allions_interfere$Isotopologue),
                                          '<br>',
                                          "Ion Formula:", CP_allions_interfere$Adduct_Formula,
                                          '<br>',
                                          "Adduct isotopes:", paste0("[12C]:", CP_allions_interfere$`12C`, "  [13C]:", CP_allions_interfere$`13C`,
                                                                     "  [35Cl]:", CP_allions_interfere$`35Cl`, "  [37Cl]:", CP_allions_interfere$`37Cl`))
                    )
                |>
                    plotly::layout(xaxis = list(title = "Number of carbons (12C+13C)"),
                                   yaxis = list(title = "Number of chlorines (35Cl+37Cl)"),
                                   legend=list(title=list(text='<b> Interference at MS res? </b>')))
            )

        }

        # Output the interference bar plot: Rel_ab vs m/z

        if ("79Br" %in% names(CP_allions_interfere) == TRUE){
            output$Plotly2 <- plotly::renderPlotly(
                p <- CP_allions_interfere |> plotly::plot_ly(
                    x = ~`m/z`,
                    y = ~Rel_ab,
                    type = "bar",
                    color = ~interference,
                    #text = ~Adduct,
                    hoverinfo = "text",
                    hovertext = paste("Molecule_Formula:", CP_allions_interfere$Molecule_Formula,
                                      '<br>',
                                      "Adduct/Fragment ion:", paste0(CP_allions_interfere$Adduct, CP_allions_interfere$Isotopologue),
                                      '<br>',
                                      "Ion Formula:", CP_allions_interfere$Adduct_Formula,
                                      '<br>',
                                      "Adduct isotopes:", paste0("[12C]:", CP_allions_interfere$`12C`, "  [13C]:", CP_allions_interfere$`13C`,
                                                                 "  [35Cl]:", CP_allions_interfere$`35Cl`, "  [37Cl]:", CP_allions_interfere$`37Cl`, " [79Br]:", CP_allions_interfere$`79Br`, " [81Br]:", CP_allions_interfere$`81Br`),
                                      '<br>',
                                      "m/z:", CP_allions_interfere$`m/z`,
                                      '<br>',
                                      "m/z diff (prev and next):", CP_allions_interfere$difflag, "&", CP_allions_interfere$difflead,
                                      '<br>',
                                      "Resolution needed (prev and next):", CP_allions_interfere$reslag, "&", CP_allions_interfere$reslead)
                )
                |>
                    plotly::layout(legend=list(title=list(text='<b> Interference at MS res? </b>')))
            )
        } else {
            output$Plotly2 <- plotly::renderPlotly(
                p <- CP_allions_interfere |> plotly::plot_ly(
                    x = ~`m/z`,
                    y = ~Rel_ab,
                    type = "bar",
                    color = ~interference,
                    #text = ~Adduct,
                    hoverinfo = "text",
                    hovertext = paste("Molecule_Formula:", CP_allions_interfere$Molecule_Formula,
                                      '<br>',
                                      "Adduct/Fragment ion:", paste0(CP_allions_interfere$Adduct, CP_allions_interfere$Isotopologue),
                                      '<br>',
                                      "Ion Formula:", CP_allions_interfere$Adduct_Formula,
                                      '<br>',
                                      "Adduct isotopes:", paste0("[12C]:", CP_allions_interfere$`12C`, "  [13C]:", CP_allions_interfere$`13C`,
                                                                 "  [35Cl]:", CP_allions_interfere$`35Cl`, "  [37Cl]:", CP_allions_interfere$`37Cl`),
                                      '<br>',
                                      "m/z:", CP_allions_interfere$`m/z`,
                                      '<br>',
                                      "m/z diff (prev and next):", CP_allions_interfere$difflag, "&", CP_allions_interfere$difflead,
                                      '<br>',
                                      "Resolution needed (prev and next):", CP_allions_interfere$reslag, "&", CP_allions_interfere$reslead)
                )
                |>
                    plotly::layout(legend=list(title=list(text='<b> Interference at MS res? </b>')))
            )

        }

        output$Table2 <- DT::renderDT({
            # Show data
            suppress_dt_size_warning(
                DT::datatable(CP_allions_interfere,
                              filter = "top", extensions = c("Buttons", "Scroller"),
                              options = list(scrollY = 650,
                                             scrollX = 500,
                                             deferRender = TRUE,
                                             scroller = TRUE,
                                             buttons = list(list(extend = "colvis", targets = 0, visible = FALSE)),
                                             dom = "lBfrtip",
                                             fixedColumns = TRUE),
                              rownames = FALSE)
            )
        }, server = TRUE)
    })
    # go2 end
    register_table_downloads(CP_allions_compl2, "download_interference_csv", "download_interference_xlsx", "CPions_interfering_ions")

    ##################################################################
    ############ go3: Skyline tab ############
    ##################################################################

    shiny::observeEvent(input$go3, {


if(input$skylineoutput == "mz"){ #Removed  skylineoutput==IonFormula since not compatible with [M-Cl]- (adduct not available in current skyline)

    if (input$QuantIon == "Most intense" & input$skyline_mode == "advanced") {
        CP_allions_skyline <- build_skyline_transition_list(
            CP_allions = CP_allions_glob_adv(),
            mode = "advanced",
            quant_ion = "Most intense",
            ms_resolution = input$MSresolution,
            strategy = input$skyline_strategy,
            preferred_qual_n = as.integer(input$skyline_qual_n)
        )
    } else if (input$QuantIon == "Most intense" & input$skyline_mode == "normal") {
        CP_allions_skyline <- build_skyline_transition_list(
            CP_allions = CP_allions_glob(),
            mode = "normal",
            quant_ion = "Most intense",
            ms_resolution = input$MSresolution,
            strategy = input$skyline_strategy,
            preferred_qual_n = as.integer(input$skyline_qual_n)
        )
    } else if (input$QuantIon == "Interference-filtered" & input$skyline_mode == "advanced") {

        req(CP_allions_compl2())
        CP_allions_skyline <- build_skyline_transition_list(
            CP_allions = CP_allions_compl2(),
            mode = "advanced",
            quant_ion = "Interference-filtered",
            ms_resolution = MSresolution(),
            strategy = input$skyline_strategy,
            preferred_qual_n = as.integer(input$skyline_qual_n)
        )
    } else if (input$QuantIon == "Interference-filtered" & input$skyline_mode == "normal") {

        req(CP_allions_compl2())
        CP_allions_skyline <- build_skyline_transition_list(
            CP_allions = CP_allions_compl2(),
            mode = "normal",
            quant_ion = "Interference-filtered",
            ms_resolution = MSresolution(),
            strategy = input$skyline_strategy,
            preferred_qual_n = as.integer(input$skyline_qual_n)
        )
    }

    CP_allions_skyline_reactive(CP_allions_skyline)
}


output$Table3 <- DT::renderDT({
    # Show data
    suppress_dt_size_warning(
        DT::datatable(CP_allions_skyline_reactive(),
                      filter = "top", extensions = c("Buttons", "Scroller"),
                      options = list(scrollY = 650,
                                     scrollX = 500,
                                     deferRender = TRUE,
                                     scroller = TRUE,
                                     buttons = list(list(extend = "colvis", targets = 0, visible = FALSE)),
                                     dom = "lBfrtip",
                                     fixedColumns = TRUE),
                      rownames = FALSE)
    ) |>
        DT::formatStyle(
            "Interference at MS Res?",
            target = "row",
            backgroundColor = DT::styleEqual(c("YES", "NO"), c("#fff3cd", NA)),
            fontWeight = DT::styleEqual(c("YES", "NO"), c("600", NA))
        )
}, server = TRUE)


})


########## go3 end

#----Outputs_End

    register_table_downloads(CP_allions_skyline_reactive, "download_iontable_csv", "download_iontable_xlsx", "Skyline_full_table")

    output[["download_skyline_csv"]] <- shiny::downloadHandler(
        filename = function() {
            "Skyline_transition_list_skyline.csv"
        },
        content = function(file) {
            data <- CP_allions_skyline_reactive()
            shiny::req(data)
            skyline_cols <- c(
                "Molecule List Name",
                "Molecule Name",
                "Precursor Charge",
                "Label Type",
                "Precursor m/z",
                "Explicit Retention Time",
                "Explicit Retention Time Window",
                "Note"
            )
            readr::write_csv(data[, intersect(skyline_cols, names(data))], file)
        }
    )

    output[["download_skyline_xlsx"]] <- shiny::downloadHandler(
        filename = function() {
            "Skyline_transition_list_skyline.xlsx"
        },
        content = function(file) {
            data <- CP_allions_skyline_reactive()
            shiny::req(data)
            skyline_cols <- c(
                "Molecule List Name",
                "Molecule Name",
                "Precursor Charge",
                "Label Type",
                "Precursor m/z",
                "Explicit Retention Time",
                "Explicit Retention Time Window",
                "Note"
            )
            openxlsx::write.xlsx(data[, intersect(skyline_cols, names(data))], file = file, overwrite = TRUE)
        }
    )

}
