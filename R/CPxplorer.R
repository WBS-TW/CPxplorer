library(shiny)
library(shinythemes)
library(DT)
library(plotly)
library(tidyverse)
library(readxl)
library(crosstalk)
library(enviPat)
library(markdown)
# library(rJava)
# library(rcdklibs)
# library(rcdk)


data("isotopes")
source("./R/getAdduct.R")
#source("./R/getSkyline.R")
source("./R/utils.R")



#--------------------------------UI function----------------------------------#

ui <- shiny::navbarPage(
        "CPXplorer",
        theme = shinythemes::shinytheme('spacelab'),
        shiny::tabPanel("Normal settings",
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        shiny::numericInput("Cmin", "C atoms min (allowed 3-40)", value = 9, min = 3, max = 40),
                                        shiny::numericInput("Cmax", "C atoms max (allowed 4-40)", value = 30, min = 4, max = 40),
                                        shiny::numericInput("Clmin", "Cl atoms min (allowed 1-15))", value = 3, min = 1, max = 15),
                                        shiny::numericInput("Clmax", "Cl atoms max (allowed 1-15)", value = 15, min = 1, max = 15),
                                        shiny::numericInput("Brmin", "Br atoms min (allowed 1-15))", value = 1, min = 1, max = 15),
                                        shiny::numericInput("Brmax", "Br atoms max (allowed 1-15)", value = 4, min = 1, max = 15),
                                        shiny::br(),
                                        selectInput("Adducts", "Add adducts/fragments",
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
                                                    selected = "[PCA-Cl]-",
                                                    multiple = TRUE,
                                                    selectize = TRUE,
                                                    width = NULL,
                                                    size = NULL),
                                        shiny::numericInput("threshold", "Isotope rel ab threshold (5-99%)", value = 5, min = 0, max = 99),
                                        shiny::actionButton("go1", "Submit", width = "100%"),
                                        width = 3),
                                shiny::mainPanel(
                                        DT::dataTableOutput("Table_norm", width = "100%")
                                )
                        )
                        )),
        shiny::tabPanel("Advanced settings",
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        shiny::numericInput("Cmin_adv", "C atoms min (allowed 3-40)", value = 9, min = 3, max = 40),
                                        shiny::numericInput("Cmax_adv", "C atoms max (allowed 4-40)", value = 30, min = 4, max = 40),
                                        shiny::numericInput("Clmin_adv", "Cl atoms min (allowed 1-15))", value = 3, min = 1, max = 15),
                                        shiny::numericInput("Clmax_adv", "Cl atoms max (allowed 1-15)", value = 15, min = 1, max = 15),
                                        shiny::numericInput("Brmin_adv", "Br atoms min (allowed 1-15))", value = 1, min = 1, max = 15),
                                        shiny::numericInput("Brmax_adv", "Br atoms max (allowed 1-15)", value = 4, min = 1, max = 15),
                                        shiny::br(),
                                        selectInput("Compclass_adv", "Compound Class",
                                                    choices = c("PCA", "PCO","BCA"),
                                                    selected = "PCA",
                                                    multiple = TRUE,
                                                    selectize = TRUE,
                                                    width = NULL,
                                                    size = NULL),
                                        selectInput("Adducts_adv", "Which adduct",
                                                    choices = c("-Cl", "-H", "-HCl", "-Cl-HCl","-Cl-2HCl", "-Cl-3HCl", "-2Cl-HCl", "+Cl","+Br"),
                                                    selected = "-Cl",
                                                    multiple = TRUE,
                                                    selectize = TRUE,
                                                    width = NULL,
                                                    size = NULL),
                                        selectInput("Charge_adv", "Which charge",
                                                    choices = c("-", "+"),
                                                    selected = "-",
                                                    multiple = FALSE,
                                                    selectize = TRUE,
                                                    width = NULL,
                                                    size = NULL),
                                        selectInput("TP_adv", "Transformation product?",
                                                    choices = c("None", "-Cl+OH", "-2Cl+2OH", "+OH", "+2OH", "+SO4"),
                                                    selected = "None",
                                                    multiple = TRUE,
                                                    selectize = TRUE,
                                                    width = NULL,
                                                    size = NULL),
                                        shiny::numericInput("threshold_adv", "Isotope rel ab threshold (0-99%)", value = 50, min = 0, max = 99),
                                        shiny::actionButton("go_adv", "Submit", width = "100%"),
                                        width = 3),
                                shiny::mainPanel(
                                        DT::dataTableOutput("Table_adv", width = "100%")
                                )
                        )
                        )),
        shiny::tabPanel("Interfering ions",
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        shiny::numericInput("MSresolution", "MS Resolution", value = 60000, min = 100, max = 5000000),
                                        shiny::radioButtons("interfere_NormAdv", label = "From Normal or Advanced settings", choices = c("normal", "advanced"), selected = "normal"),
                                        shiny::actionButton("go2", "Calculate", width = "100%"),
                                        width = 2
                                ),
                                shiny::mainPanel(
                                        plotly::plotlyOutput("Plotly"),
                                        plotly::plotlyOutput("Plotly2"),
                                        DT::dataTableOutput("Table2", width = "100%")
                                        
                                )
                        )
                        )),
        shiny::tabPanel("Skyline",
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        #shiny::numericInput("MSresolution2", "MS Resolution", value = 60000, min = 100, max = 3000000),
                                        shiny::radioButtons("QuantIon", label = "Use as Quant Ion", choices = c("Most intense")),
                                        #shiny::radioButtons("skylineoutput", label = "Output table", choices = c("mz", "IonFormula")),
                                        shiny::radioButtons("skylineoutput", label = "Output table", choices = c("mz")),
                                        shiny::radioButtons("skyline_NormAdv", label = "From Normal or Advanced settings", choices = c("normal", "advanced"), selected = "normal"),
                                        shiny::actionButton("go3", "Transition List", width = "100%"),
                                        width = 4
                                ),
                                shiny::mainPanel(
                                        DT::dataTableOutput("Table3", width = "100%")
                                        
                                )
                        )
                        )),
        shiny::tabPanel(
                "Instructions",
                shiny::sidebarLayout(
                        shiny::sidebarPanel(shiny::h3("Manual"),
                                            width = 3),
                        shiny::mainPanel(
                                shiny::includeMarkdown("instructions.md")
                        )
                )
        )
)



#---------------------------Shiny Server function------------------------------#



server = function(input, output, session) {
        
        # Set reactive values from user input
        
        # GENERAL
        MSresolution <- eventReactive(input$go2, {as.integer(input$MSresolution)})
        
        # NORMAL
        selectedAdducts <- eventReactive(input$go1, {as.character(input$Adducts)})
        
        C <- eventReactive(input$go1, {as.integer(input$Cmin:input$Cmax)})
        Cl <- eventReactive(input$go1, {as.integer(input$Clmin:input$Clmax)})
        Clmin <- eventReactive(input$go1, {as.integer(input$Clmin)})
        Clmax <- eventReactive(input$go1, {as.integer(input$Clmax)})
        Br <- eventReactive(input$go1, {as.integer(input$Brmin:input$Brmax)})
        Brmin <- eventReactive(input$go1, {as.integer(input$Brmin)})
        Brmax <- eventReactive(input$go1, {as.integer(input$Brmax)})
        threshold <- eventReactive(input$go1, {as.integer(input$threshold)})
        
        
        # ADVANCED
        selectedCompounds_adv <- shiny::eventReactive(input$go_adv, {as.character((input$Compclass_adv))})
        selectedAdducts_adv <- shiny::eventReactive(input$go_adv, {as.character(input$Adducts_adv)})
        selectedCharge_adv <- shiny::eventReactive(input$go_adv, {as.character((input$Charge_adv))})
        selectedTP_adv <- shiny::eventReactive(input$go_adv, {as.character((input$TP_adv))})
        
        C_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Cmin:input$Cmax_adv)})
        Cl_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Clmin:input$Clmax_adv)})
        Clmin_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Clmin_adv)})
        Clmax_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Clmax_adv)})
        Br_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Brmin:input$Brmax_adv)})
        Brmin_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Brmin_adv)})
        Brmax_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$Brmax_adv)})
        threshold_adv <- shiny::eventReactive(input$go_adv, {as.integer(input$threshold_adv)})
        
        
        #----Outputs_Start
        # NORMAL
        CP_allions_glob <- eventReactive(input$go1, {
                
                # Create a Progress bar object
                progress <- shiny::Progress$new()
                
                # Make sure it closes when we exit this reactive, even if there's an error
                on.exit(progress$close())
                progress$set(message = "Calculating", value = 0)
                
                Adducts <- as.character(selectedAdducts())
                
                # function to get adducts or fragments
                CP_allions <- data.frame(Molecule_Formula = character(), Halo_perc = double())
                for (i in seq_along(Adducts)) {
                        progress$inc(1/length(Adducts), detail = paste0("Adduct: ", Adducts[i], " . Please wait.."))
                        if(str_detect(Adducts[i], "\\bBCA\\b")){
                                input <- getAdduct_BCA(adduct_ions = Adducts[i], C = C(), Cl = Cl(), Clmax = Clmax(), 
                                                       Br = Br(), Brmax = Brmax(), threshold = threshold())
                        } else {
                                input <- getAdduct_normal(adduct_ions = Adducts[i], C = C(), Cl = Cl(), Clmax = Clmax(), threshold = threshold())
                        }
                        CP_allions <- dplyr::full_join(CP_allions, input)
                }
                return(CP_allions)
        })
        
        #############################################################
        ### go1: Calculate the isotopes from initial settings tab ###
        #############################################################
        
        shiny::observeEvent(input$go1, {
                output$Table_norm <- DT::renderDT(server=FALSE,{ #need to keep server = FALSE otherwise excel download the visible rows of the table, but this will also give warning about large tables
                        # Show data
                        DT::datatable(CP_allions_glob(), 
                                      filter = "top", extensions = c("Buttons", "Scroller"),
                                      options = list(scrollY = 650,
                                                     scrollX = 500,
                                                     deferRender = TRUE,
                                                     scroller = TRUE,
                                                     buttons = list(list(extend = "excel", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "csv", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "colvis", targets = 0, visible = FALSE)),
                                                     dom = "lBfrtip",
                                                     fixedColumns = TRUE), 
                                      rownames = FALSE)
                })
        })
        # go1 end
        
        
        # ADVANCED
        CP_allions_glob_adv <- eventReactive(input$go_adv, {
                
                # Create a Progress bar object
                progress <- shiny::Progress$new()
                
                # Make sure it closes when we exit this reactive, even if there's an error
                on.exit(progress$close())
                progress$set(message = "Calculating", value = 0)
                
                Compounds <- as.character(selectedCompounds_adv())
                Adducts <- as.character(selectedAdducts_adv())
                Charge <- as.character(selectedCharge_adv())
                TP <- as.character(selectedTP_adv())
                
                # function to get adducts or fragments
                CP_allions <- data.frame(Molecule_Formula = character(), Molecule_Halo_perc = double())
                
                # nested for loop to get all combinations of Compounds, Adducts, TP
                for (i in seq_along(Compounds)) {
                        progress$inc(1/length(Compounds), detail = paste0("Compound Class: ", Compounds[i], " . Please wait.."))
                        
                        for (j in seq_along(Adducts)) {
                                
                                for (k in seq_along(TP)) {
                                        input <- getAdduct_advanced(Compounds = Compounds[i], Adduct_Ion = Adducts[j], TP = TP[k], Charge = Charge,
                                                                    C = C_adv(), Cl = Cl_adv(), Clmax = Clmax_adv(), Br = Br_adv(), Brmax = Brmax_adv(), 
                                                                    threshold = threshold_adv())
                                        CP_allions <- dplyr::full_join(CP_allions, input)
                                }
                        }
                }
                return(CP_allions)
        })
        
        ### go_adv: Calculate the isotopes from initial settings tab ###
        shiny::observeEvent(input$go_adv, {
                output$Table_adv <- DT::renderDT(server=FALSE,{ #need to keep server = FALSE otherwise excel download the visible rows of the table, but this will also give warning about large tables
                        # Show data
                        DT::datatable(CP_allions_glob_adv(), 
                                      filter = "top", extensions = c("Buttons", "Scroller"),
                                      options = list(scrollY = 650,
                                                     scrollX = 500,
                                                     deferRender = TRUE,
                                                     scroller = TRUE,
                                                     buttons = list(list(extend = "excel", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "csv", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "colvis", targets = 0, visible = FALSE)),
                                                     dom = "lBfrtip",
                                                     fixedColumns = TRUE), 
                                      rownames = FALSE)
                })
        })
        # go_adv end
        
        ##################################################################         
        ############ go2: Calculates the interfering ions tab ############
        ##################################################################         
        
        shiny::observeEvent(input$go2, {
                
                if (input$interfere_NormAdv == "normal") {
                        CP_allions_compl2 <- CP_allions_glob()}
                else if (input$interfere_NormAdv == "advanced") {CP_allions_compl2 <- CP_allions_glob_adv()}
                
                CP_allions_compl2 <- CP_allions_compl2 |>
                        arrange(`m/z`) |>
                        mutate(difflag = round(abs(`m/z` - lag(`m/z`, default = first(`m/z`))),6)) |>
                        mutate(difflead = round(abs(`m/z` - lead(`m/z`, default = last(`m/z`))), 6)) |>
                        mutate(reslag = round(`m/z`/difflag, 0)) |>
                        mutate(reslead = round(`m/z`/difflead, 0)) |>
                        mutate(interference = case_when(
                                difflag == 0 | difflead == 0 ~ TRUE, # need to keep this true to make same mass ions TRUE
                                reslag >= as.integer(MSresolution()) | reslead >= as.integer(MSresolution()) ~ TRUE,
                                reslag < as.integer(MSresolution()) & reslead < as.integer(MSresolution()) ~ FALSE
                        )
                        )
                # change first and last row to false since their lead/lag is zero
                CP_allions_compl2$interference[1] <- FALSE
                CP_allions_compl2$interference[length(CP_allions_compl2$interference)] <- FALSE
                
                # Output scatterplot: #Cl vs #C  if Br exists
                if ("79Br" %in% names(CP_allions_compl2) == TRUE){ 
                        output$Plotly <- plotly::renderPlotly(
                                p <- CP_allions_compl2 |> 
                                        mutate(`79Br` = tidyr::replace_na(`79Br`, 0)) |>
                                        mutate(`81Br` = tidyr::replace_na(`81Br`, 0)) |>
                                        plot_ly(
                                                x = ~ (`12C`+`13C`), 
                                                y = ~(`35Cl`+`37Cl`+`79Br`+`81Br`), # need to incorporate Br later
                                                type = "scatter",
                                                mode = "markers",
                                                color = ~interference,
                                                hoverinfo = "text",
                                                hovertext = paste("Molecule_Formula:", CP_allions_compl2$Molecule_Formula,
                                                                  '<br>',
                                                                  "Adduct/Fragment ion:", CP_allions_compl2$Adduct,
                                                                  '<br>',
                                                                  "Ion Formula:", CP_allions_compl2$Adduct_Formula,
                                                                  '<br>',
                                                                  "Adduct isotopes:", paste0("[12C]:", CP_allions_compl2$`12C`, "  [13C]:", CP_allions_compl2$`13C`, 
                                                                                             "  [35Cl]:", CP_allions_compl2$`35Cl`, "  [37Cl]:", CP_allions_compl2$`37Cl`, " [79Br]:", CP_allions_compl2$`79Br`, " [81Br]:", CP_allions_compl2$`81Br`))
                                        )
                                |> 
                                        plotly::layout(xaxis = list(title = "Number of carbons (12C+13C)"),
                                                       yaxis = list(title = "Number of halogens (35Cl+37Cl+79Br+81Br)"),
                                                       legend=list(title=list(text='<b> Interference at MS res? </b>')))
                        )
                } else { #if there are no bromines
                        output$Plotly <- plotly::renderPlotly(
                                p <- CP_allions_compl2 |> 
                                        plot_ly(
                                                x = ~ (`12C`+`13C`), 
                                                y = ~(`35Cl`+`37Cl`),
                                                type = "scatter",
                                                mode = "markers",
                                                color = ~interference,
                                                hoverinfo = "text",
                                                hovertext = paste("Molecule_Formula:", CP_allions_compl2$Molecule_Formula,
                                                                  '<br>',
                                                                  "Adduct/Fragment ion:", CP_allions_compl2$Adduct,
                                                                  '<br>',
                                                                  "Ion Formula:", CP_allions_compl2$Adduct_Formula,
                                                                  '<br>',
                                                                  "Adduct isotopes:", paste0("[12C]:", CP_allions_compl2$`12C`, "  [13C]:", CP_allions_compl2$`13C`, 
                                                                                             "  [35Cl]:", CP_allions_compl2$`35Cl`, "  [37Cl]:", CP_allions_compl2$`37Cl`))
                                        )
                                |> 
                                        plotly::layout(xaxis = list(title = "Number of carbons (12C+13C)"),
                                                       yaxis = list(title = "Number of chlorines (35Cl+37Cl)"),
                                                       legend=list(title=list(text='<b> Interference at MS res? </b>')))
                        )
                        
                }
                
                # Output the interference bar plot: Rel_ab vs m/z
                
                if ("79Br" %in% names(CP_allions_compl2) == TRUE){ 
                        output$Plotly2 <- plotly::renderPlotly(
                                p <- CP_allions_compl2 |> plot_ly(
                                        x = ~`m/z`, 
                                        y = ~Rel_ab,
                                        type = "bar",
                                        color = ~interference,
                                        #text = ~Adduct,
                                        hoverinfo = "text",
                                        hovertext = paste("Molecule_Formula:", CP_allions_compl2$Molecule_Formula,
                                                          '<br>',
                                                          "Adduct/Fragment ion:", CP_allions_compl2$Adduct,
                                                          '<br>',
                                                          "Ion Formula:", CP_allions_compl2$Adduct_Formula,
                                                          '<br>',
                                                          "Adduct isotopes:", paste0("[12C]:", CP_allions_compl2$`12C`, "  [13C]:", CP_allions_compl2$`13C`, 
                                                                                     "  [35Cl]:", CP_allions_compl2$`35Cl`, "  [37Cl]:", CP_allions_compl2$`37Cl`, " [79Br]:", CP_allions_compl2$`79Br`, " [81Br]:", CP_allions_compl2$`81Br`),
                                                          '<br>',
                                                          "m/z:", CP_allions_compl2$`m/z`,
                                                          '<br>',
                                                          "m/z diff (prev and next):", CP_allions_compl2$difflag, "&", CP_allions_compl2$difflead,
                                                          '<br>',
                                                          "Resolution needed (prev and next):", CP_allions_compl2$reslag, "&", CP_allions_compl2$reslead)
                                )
                                |> 
                                        plotly::layout(legend=list(title=list(text='<b> Interference at MS res? </b>')))
                        )
                } else {
                        output$Plotly2 <- plotly::renderPlotly(
                                p <- CP_allions_compl2 |> plot_ly(
                                        x = ~`m/z`, 
                                        y = ~Rel_ab,
                                        type = "bar",
                                        color = ~interference,
                                        #text = ~Adduct,
                                        hoverinfo = "text",
                                        hovertext = paste("Molecule_Formula:", CP_allions_compl2$Molecule_Formula,
                                                          '<br>',
                                                          "Adduct/Fragment ion:", CP_allions_compl2$Adduct,
                                                          '<br>',
                                                          "Ion Formula:", CP_allions_compl2$Adduct_Formula,
                                                          '<br>',
                                                          "Adduct isotopes:", paste0("[12C]:", CP_allions_compl2$`12C`, "  [13C]:", CP_allions_compl2$`13C`, 
                                                                                     "  [35Cl]:", CP_allions_compl2$`35Cl`, "  [37Cl]:", CP_allions_compl2$`37Cl`),
                                                          '<br>',
                                                          "m/z:", CP_allions_compl2$`m/z`,
                                                          '<br>',
                                                          "m/z diff (prev and next):", CP_allions_compl2$difflag, "&", CP_allions_compl2$difflead,
                                                          '<br>',
                                                          "Resolution needed (prev and next):", CP_allions_compl2$reslag, "&", CP_allions_compl2$reslead)
                                )
                                |> 
                                        plotly::layout(legend=list(title=list(text='<b> Interference at MS res? </b>')))
                        )
                        
                }
                
                output$Table2 <- DT::renderDT(server=FALSE,{ #need to keep server = FALSE otherwise excel download only part of rows
                        # Show data
                        DT::datatable(CP_allions_compl2, 
                                      filter = "top", extensions = c("Buttons", "Scroller"),
                                      options = list(scrollY = 650,
                                                     scrollX = 500,
                                                     deferRender = TRUE,
                                                     scroller = TRUE,
                                                     buttons = list(list(extend = "excel", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "csv", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "colvis", targets = 0, visible = FALSE)),
                                                     dom = "lBfrtip",
                                                     fixedColumns = TRUE), 
                                      rownames = FALSE)
                })
        })
        # go2 end
        
        ##################################################################         
        ############ go3: Skyline tab ############
        ################################################################## 
        
        shiny::observeEvent(input$go3, {
                
                
                if(input$skylineoutput == "mz"){ #Removed  skylineoutput==IonFormula since not compatible with [M-Cl]- (adduct not available in skyline)
                        
                        if (input$skyline_NormAdv == "advanced") {
                                CP_allions_skyline <- CP_allions_glob_adv() |>
                                        mutate(`Molecule List Name` = case_when(Compound_Class == "PCA" ~ paste0("PCA-C", str_extract(Molecule_Formula, "(?<=C)\\d+(?=H)")),
                                                                                Compound_Class == "PCO" ~ paste0("PCO-C", str_extract(Molecule_Formula, "(?<=C)\\d+(?=H)")),
                                                                                Compound_Class == "BCA" ~ paste0("BCA-C", str_extract(Molecule_Formula, "(?<=C)\\d+(?=H)")))) |>
                                        rename(`Molecule Name` = Molecule_Formula) |>
                                        mutate(`Precursor m/z` = `m/z`) |>
                                        # mutate(Note = str_replace(Adduct, "\\].*", "]")) |>
                                        # mutate(Note = str_replace(Note, "(.+?(?=\\-))|(.+?(?=\\+))", "[M")) |>
                                        mutate(Note = ifelse(TP == "None", Compound_Class, paste0(Compound_Class, TP))) |>
                                        rename(`Precursor Charge` = Charge) |>
                                        add_column(`Explicit Retention Time` = NA) |>
                                        add_column(`Explicit Retention Time Window` = NA) |>
                                        group_by(`Molecule Name`) |>
                                        mutate(`Label Type` = ifelse(Rel_ab == 100, "Quan", "Qual")) |> # choose the highest rel_ab ion as quan ion and the rest will be qual
                                        ungroup() |>
                                        select(`Molecule List Name`,
                                               `Molecule Name`,
                                               `Precursor Charge`,
                                               `Label Type`,
                                               `Precursor m/z` = `m/z`,
                                               `Explicit Retention Time`,
                                               `Explicit Retention Time Window`,
                                               Note)
                        } else if (input$skyline_NormAdv == "normal") {
                                CP_allions_skyline <- CP_allions_glob() |> 
                                        mutate(`Molecule List Name` = case_when(str_detect(Adduct, "(?<=.)PCA(?=.)") == TRUE ~ paste0("PCA-C", str_extract(Molecule_Formula, "(?<=C)\\d+(?=H)")),
                                                                                str_detect(Adduct, "(?<=.)PCO(?=.)") == TRUE ~ paste0("PCO-C", str_extract(Molecule_Formula, "(?<=C)\\d+(?=H)")),
                                                                                str_detect(Adduct, "(?<=.)BCA(?=.)") == TRUE ~ paste0("BCA-C", str_extract(Molecule_Formula, "(?<=C)\\d+(?=H)")))) |>
                                        rename(`Molecule Name` = Molecule_Formula) |>
                                        mutate(`Precursor m/z` = `m/z`) |>
                                        # mutate(Note = str_replace(Adduct, "\\].*", "]")) |>
                                        # mutate(Note = str_replace(Note, "(.+?(?=\\-))|(.+?(?=\\+))", "[M")) |>
                                        rename(Note = Adduct) |>
                                        rename(`Precursor Charge` = Charge) |>
                                        add_column(`Explicit Retention Time` = NA) |>
                                        add_column(`Explicit Retention Time Window` = NA) |>
                                        group_by(`Molecule Name`) |>
                                        mutate(`Label Type` = ifelse(Rel_ab == 100, "Quan", "Qual")) |> # choose the highest rel_ab ion as quan ion and the rest will be qual
                                        ungroup() |>
                                        select(`Molecule List Name`,
                                               `Molecule Name`,
                                               `Precursor Charge`,
                                               `Label Type`,
                                               `Precursor m/z` = `m/z`,
                                               `Explicit Retention Time`,
                                               `Explicit Retention Time Window`,
                                               Note)
                        }
                }
                
                
                output$Table3 <- DT::renderDT(server=FALSE,{ #need to keep server = FALSE otherwise excel download the visible rows of the table, but this will also give warning about large tables
                        # Show data
                        DT::datatable(CP_allions_skyline,
                                      filter = "top", extensions = c("Buttons", "Scroller"),
                                      options = list(scrollY = 650,
                                                     scrollX = 500,
                                                     deferRender = TRUE,
                                                     scroller = TRUE,
                                                     buttons = list(list(extend = "excel", filename = "Skyline_transition_list", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "csv", filename = "Skyline_transition_list", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "colvis", targets = 0, visible = FALSE)),
                                                     dom = "lBfrtip",
                                                     fixedColumns = TRUE),
                                      rownames = FALSE)
                })
                
                
        })
        
        
        ########## go3 end
        
        #----Outputs_End
        
        
        
        
        
        
        # Close the app when the session ends
        if(!interactive()) {
                session$onSessionEnded(function() {
                        stopApp()
                        q("no")
                })
        }
        
}

shiny::shinyApp(ui, server)
