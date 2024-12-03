# CPquant: Shiny CPs Quantification for Skyline Output
#
#
#
library(shiny)
library(tidyverse)
library(readxl)
library(plotly)
library(DT)
library(nnls)
library(openxlsx)



options(shiny.maxRequestSize = 15000 * 1024^2)


# UI
ui <- shiny::navbarPage("Quantification by deconvolution from Skyline output",
                        shiny::tabPanel("Quantification Inputs", 
                                        shiny::fluidPage(shiny::sidebarLayout(
                                                shiny::sidebarPanel(
                                                        shiny::fileInput("fileInput", "Import excel file from Skyline", 
                                                                         accept = c('xlsx')),
                                                        # shiny::radioButtons("quantStd",
                                                        #                    label = "Quantification based on chain length or mixture standards?",
                                                        #                   choices = c("Chain length",
                                                        #                              "Mixture",
                                                        #                             "Both")),
                                                        shiny::radioButtons("blankSubtraction",
                                                                            label = "Blank signal subtraction",
                                                                            choices = c("Yes", "No"), selected = "No") ,
                                                        # shiny::selectInput("includedCPs", "Include which CPs for quantification?",
                                                        #                   choices = c("vSCCPs", "SCCPs", "MCCPs", "LCCPs", "vLCCPs"),
                                                        #                  selected = c("vSCCPs", "SCCPs", "MCCPs", "LCCPs", "vLCCPs"),
                                                        #                 multiple = TRUE),
                                                        # New radio button to select "IS only"
                                                        shiny::radioButtons("includeISOnly", label = "Correct the signal with RS?",
                                                                            choices = c("Yes", "No"), selected = "No") ,
                                                        
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
                                        DT::DTOutput("LOD")        # Second table output (LOD table)
                                )
                        )
                        
)


################################################################################
server <- function(input, output, session) {
        
        # Create a function to process each row
        
        Skyline_output <- reactive({
                req(input$fileInput) #requires that the input is available
                df <- readxl::read_excel(input$fileInput$datapath) |> 
                        mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) |> 
                        mutate(Area = as.numeric(Area)) |> 
                        mutate(Area = replace_na(Area, 0)) |>  # Replace missing values with 0
                        #mutate(RatioQuanToQual = as.numeric(RatioQuanToQual)) |> 
                        #mutate(RatioQualToQuan = as.numeric(RatioQualToQuan)) |> 
                        # Extract significant parts from the 'Molecule' column
                        # Extract "C" followed by numbers
                        mutate(C_part = str_extract(Molecule, "C\\d+"),
                               
                               # Extract "Cl" followed by numbers
                               Cl_part = str_extract(Molecule, "Cl\\d+"),
                               
                               # Extract numeric values for sorting
                               C_number = as.numeric(str_extract(C_part, "\\d+")),
                               Cl_number = as.numeric(str_extract(Cl_part, "\\d+")),
                               
                               # Combine them into SimplifiedMolecule
                               SimplifiedMolecule = str_c(C_part, Cl_part, sep = "")
                        ) 
                
                
                
                #  Normalize data based on 'Include IS only' input using case_when
                df <- df |> 
                        group_by(`Replicate Name`) |>
                        filter(`Isotope Label Type` == "Quan") |>
                        mutate(Area = case_when(
                                input$includeISOnly == "Yes" & any(Molecule == "RS") ~ Area / first(Area[Molecule == "RS"]),
                                TRUE ~ Area
                        )) |> 
                        ungroup()
                
                df
                
                
                
                # Calculate the average blank value
                dfb <- df |>
                        filter(`Sample Type` == "Blank") |>
                        filter(`Isotope Label Type` == "Quan") |>
                        summarize(AverageBlank = mean(Area, na.rm = TRUE)) %>%
                        pull(AverageBlank)  # Convert to a single value
                
                # Create a logical vector for the condition
                apply_blank_subtraction <- input$blankSubtraction == "Yes"
                
                # Mutate the DataFrame with blank subtraction and handle negative values
                df <- df |>
                        mutate(
                                Area = pmax(0, Area - if (apply_blank_subtraction) dfb else 0)
                        ) |>
                        ungroup()
                
                
                
        })
        
        
        
        ###START: Define input variables       
        output$defineVariables <- shiny::renderUI({
                data <- Skyline_output()
                
                
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
                                        selected = "Note"
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
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::selectInput(
                                        inputId = "removeSamples", #select if some samples will be removed from quantification 
                                        label = 'Samples to remove from quantification?',
                                        choices = unique(Skyline_output()$`Replicate Name`),
                                        multiple = TRUE
                                )
                        ),
                        # shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), 
                        #  shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        # shiny::column(
                        # 6,
                        #shiny::sliderInput(
                        #       inputId = "removeAreas", #remove low peak areas
                        #      label = "Keep absolute peak areas above this threshold (0 means keep everything)",
                        #     min = min(Skyline_output()$Area),
                        #    max = max(Skyline_output()$Area),
                        #   value = 0,
                        #  step = 100
                        # )
                        #),
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::sliderInput(
                                        inputId = "removeRsquared", #keep only Molecule above this rsquared, zero means keep everything 
                                        label = 'Keep the the calibration curves that show rsquared above this threshold (0 means keep everything)',
                                        min = 0,
                                        max = 1,
                                        value = 0,
                                        step = 0.05
                                )
                        ),
                        #shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), 
                        #shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        #shiny::column(
                        #       6,
                        #      shiny::radioButtons("ISRS", label = "Included IS and RS?",
                        #                         choices = c("None", "IS only", "RS only", "Both IS and RS"))
                        #  )      
                )
        })
        ### END: Define input variables
        
        # Set reactive values from user input
        # removeAreas <- eventReactive(input$go, {as.numeric(input$removeAreas)})
        removeRsquared <- eventReactive(input$go, {as.numeric(input$removeRsquared)})
        
        
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
                        # Extract numeric parts for C and Cl
                        dplyr::mutate(
                                C_number = as.numeric(stringr::str_extract(Molecule, "(?<=C)\\d+")),
                                Cl_number = as.numeric(stringr::str_extract(Molecule, "(?<=Cl)\\d+"))
                        ) |>
                        # Create a composite ordering factor
                        dplyr::mutate(
                                OrderedMolecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)]))
                        )
                
                plotly::plot_ly(
                        data,
                        x = ~ OrderedMolecule,
                        y = ~ Area,
                        color = ~ `Sample Type`,
                        type = "box"
                )
        })
        
        #output$CalibrationSCCPs <- plotly::renderPlotly({
        # data <- Skyline_output() |>
        #  dplyr::filter(`Isotope Label Type` == "Quan") |> 
        # dplyr::filter(stringr::str_detect(Note, "S-")) |> 
        #dplyr::filter(`Molecule List` %in% c("PCA-C10", "PCA-C11", "PCA-C12", "PCA-C13")) |> 
        #dplyr::group_by(Note) |> 
        #dplyr::group_by(Molecule) 
        
        #  plotly::plot_ly(data, x = ~ `Analyte Concentration`, y = ~ Area, z = ~ Note, color = ~ Molecule, type = "scatter", mode = "lines+markers", 
        #                 text = ~ paste("Homologue:", Molecule, "<br>Area:", Area, "<br>Analyte Concentration (ug/g):", `Analyte Concentration`, "<br>Standard:", `Note`),
        #                hoverinfo = "text") |>
        # plotly::layout(
        #  title = "Calibration PCAs-C10-13",
        # xaxis = list(title = "Analyte Concentration"),
        #yaxis = list(title = "Area")
        #)
        #})
        
        output$CalibrationSCCPs <- plotly::renderPlotly({
                data <- Skyline_output() |>
                        dplyr::filter(`Isotope Label Type` == "Quan") |>
                        dplyr::filter(stringr::str_detect(Note, "S-")) |>
                        dplyr::filter(`Molecule List` %in% c("PCA-C10", "PCA-C11", "PCA-C12", "PCA-C13")) |>
                        dplyr::group_by(Note, `Molecule`)
                
                # Fit linear models for each combination of Note and Molecule List
                models <- data |>
                        dplyr::group_modify(~ {
                                lm_model <- lm(Area ~ `Analyte Concentration`, data = .x)
                                .x$predicted <- predict(lm_model, newdata = .x)
                                .x
                        })
                
                # Initialize the plot
                p <- plotly::plot_ly()
                
                # Add scatter plot markers with visibility set to "legendonly"
                p <- p |>
                        plotly::add_markers(
                                data = models,
                                x = ~ `Analyte Concentration`,
                                y = ~ Area,
                                color = ~ `Molecule`,
                                symbol = ~ Note, # Differentiates by Note using different symbols
                                text = ~ paste("Homologue:", `Molecule`, "<br>Area:", Area,
                                               "<br>Analyte Concentration (ug/g):", `Analyte Concentration`,
                                               "<br>Standard:", Note),
                                hoverinfo = "text",
                                visible = "legendonly",  # Initially hide this trace
                                legendgroup = ~ `Molecule` # Group by Molecule List
                        )
                
                # Add linear model lines with visibility set to "legendonly"
                p <- p |>
                        plotly::add_lines(
                                data = models,
                                x = ~ `Analyte Concentration`,
                                y = ~ predicted,
                                color = ~ `Molecule`,
                                line = list(dash = "solid"),
                                hoverinfo = "none",
                                visible = "legendonly",  # Initially hide this trace
                                legendgroup = ~ `Molecule` # Group by Molecule List
                        )
                
                # Layout configuration
                p |>
                        plotly::layout(
                                title = "Calibration PCAs-C10-13",
                                xaxis = list(title = "Analyte Concentration"),
                                yaxis = list(title = "Area")
                        )
        })
        
        
        
        
        # output$CalibrationMCCPs <- plotly::renderPlotly({
        #  data <- Skyline_output() |>
        #   dplyr::filter(`Isotope Label Type` == "Quan") |> 
        #  dplyr::filter(stringr::str_detect(Note, "M-")) |> 
        # dplyr::filter(`Molecule List` %in% c("PCA-C14", "PCA-C15", "PCA-C16", "PCA-C17")) |> 
        #  dplyr::group_by(Note) |> 
        # dplyr::group_by(Molecule)
        
        #  plotly::plot_ly(data, x = ~ `Analyte Concentration`, y = ~ Area, color = ~ Molecule, type = "scatter", mode = "lines+markers", 
        #                 text = ~ paste("Homologue:", Molecule, "<br>Area:", Area, "<br>Analyte Concentration (ug/g):", `Analyte Concentration`, "<br>Standard:", `Note`),
        #                hoverinfo = "text") |>
        # plotly::layout(
        #  title = "Calibration PCAs-C14-17",
        # xaxis = list(title = "Analyte Concentration"),
        #yaxis = list(title = "Area")
        #)
        #})
        output$CalibrationMCCPs <- plotly::renderPlotly({
                # Filter and group the data
                data <- Skyline_output() |>
                        dplyr::filter(`Isotope Label Type` == "Quan") |> 
                        dplyr::filter(stringr::str_detect(Note, "M-")) |> 
                        dplyr::filter(`Molecule List` %in% c("PCA-C14", "PCA-C15", "PCA-C16", "PCA-C17")) |> 
                        dplyr::group_by(Note, `Molecule`)
                
                # Fit linear models for each combination of Note and Molecule
                models <- data |>
                        dplyr::group_modify(~ {
                                lm_model <- lm(Area ~ `Analyte Concentration`, data = .x)
                                .x$predicted <- predict(lm_model, newdata = .x)
                                .x
                        })
                
                # Initialize the plot
                p <- plotly::plot_ly()
                
                # Add scatter plot markers with visibility set to "legendonly"
                p <- p |>
                        plotly::add_markers(
                                data = models,
                                x = ~ `Analyte Concentration`,
                                y = ~ Area,
                                color = ~ `Molecule`,
                                symbol = ~ Note, # Differentiates by Note using different symbols
                                text = ~ paste("Homologue:", `Molecule`, "<br>Area:", Area,
                                               "<br>Analyte Concentration (ug/g):", `Analyte Concentration`,
                                               "<br>Standard:", Note),
                                hoverinfo = "text",
                                visible = "legendonly",  # Initially hide this trace
                                legendgroup = ~ `Molecule` # Group by Molecule
                        )
                
                # Add linear model lines with visibility set to "legendonly"
                p <- p |>
                        plotly::add_lines(
                                data = models,
                                x = ~ `Analyte Concentration`,
                                y = ~ predicted,
                                color = ~ `Molecule`,
                                line = list(dash = "solid"),
                                hoverinfo = "none",
                                visible = "legendonly",  # Initially hide this trace
                                legendgroup = ~ `Molecule` # Group by Molecule
                        )
                
                # Layout configuration
                p |>
                        plotly::layout(
                                title = "Calibration PCAs-C14-17",
                                xaxis = list(title = "Analyte Concentration"),
                                yaxis = list(title = "Area")
                        )
        })
        
        
        #output$CalibrationLCCPs <- plotly::renderPlotly({
        # data <- Skyline_output() |>
        #  dplyr::filter(`Isotope Label Type` == "Quan") |> 
        # dplyr::filter(stringr::str_detect(Note, "L-")) |> 
        #dplyr::filter(`Molecule List` %in% c("PCA-C18", "PCA-C19", "PCA-C20", "PCA-C21", "PCA-C22", "PCA-C23", "PCA-C24", "PCA-C25", "PCA-C26", "PCA-C27", "PCA-C28", "PCA-C29", "PCA-C30")) |> 
        #dplyr::group_by(Note) |> 
        #dplyr::group_by(Molecule)
        
        #plotly::plot_ly(data, x = ~ `Analyte Concentration`, y = ~ Area, color = ~ Molecule, type = "scatter", mode = "lines+markers", 
        #               text = ~ paste("Homologue:", Molecule, "<br>Area:", Area, "<br>Analyte Concentration (ug/g):", `Analyte Concentration`, "<br>Standard:", `Note`),
        #              hoverinfo = "text") |>
        #plotly::layout(
        # title = "Calibration PCAs-C18-30",
        #xaxis = list(title = "Analyte Concentration"),
        #yaxis = list(title = "Area")
        #)
        #})
        
        output$CalibrationLCCPs <- plotly::renderPlotly({
                # Filter and group the data
                data <- Skyline_output() |>
                        dplyr::filter(`Isotope Label Type` == "Quan") |> 
                        dplyr::filter(stringr::str_detect(Note, "L-")) |> 
                        dplyr::filter(`Molecule List` %in% c("PCA-C18", "PCA-C19", "PCA-C20", "PCA-C21", "PCA-C22", "PCA-C23", 
                                                             "PCA-C24", "PCA-C25", "PCA-C26", "PCA-C27", "PCA-C28", "PCA-C29", "PCA-C30")) |> 
                        dplyr::group_by(Note, `Molecule`)
                
                # Fit linear models for each combination of Note and Molecule
                models <- data |>
                        dplyr::group_modify(~ {
                                lm_model <- lm(Area ~ `Analyte Concentration`, data = .x)
                                .x$predicted <- predict(lm_model, newdata = .x)
                                .x
                        })
                
                # Initialize the plot
                p <- plotly::plot_ly()
                
                # Add scatter plot markers with visibility set to "legendonly"
                p <- p |>
                        plotly::add_markers(
                                data = models,
                                x = ~ `Analyte Concentration`,
                                y = ~ Area,
                                color = ~ `Molecule`,
                                symbol = ~ Note, # Differentiates by Note using different symbols
                                text = ~ paste("Homologue:", `Molecule`, "<br>Area:", Area,
                                               "<br>Analyte Concentration (ug/g):", `Analyte Concentration`,
                                               "<br>Standard:", Note),
                                hoverinfo = "text",
                                visible = "legendonly",  # Initially hide this trace
                                legendgroup = ~ `Molecule` # Group by Molecule
                        )
                
                # Add linear model lines with visibility set to "legendonly"
                p <- p |>
                        plotly::add_lines(
                                data = models,
                                x = ~ `Analyte Concentration`,
                                y = ~ predicted,
                                color = ~ `Molecule`,
                                line = list(dash = "solid"),
                                hoverinfo = "none",
                                visible = "legendonly",  # Initially hide this trace
                                legendgroup = ~ `Molecule` # Group by Molecule
                        )
                
                # Layout configuration
                p |>
                        plotly::layout(
                                title = "Calibration PCAs-C18-30",
                                xaxis = list(title = "Analyte Concentration"),
                                yaxis = list(title = "Area")
                        )
        })
        
        
        
        ##### START: Deconvolution script
        
        shiny::observeEvent(input$go, {
                
                Skyline_output_filt <- Skyline_output() #|> 
                #filter(Area >= removeAreas())
                #mutate(Area = if_else(Area <= removeAreas(), 0, Area))
                
                #browser()
                CPs_standards <- Skyline_output_filt |> 
                        filter(`Sample Type` == "Standard",
                               Molecule != "IS",
                               Molecule != "RS",
                               `Isotope Label Type` == "Quan",
                               Note != "NA") |> 
                        group_by(!!input$standardAnnoColumn, Molecule) |>
                        mutate(rel_int = Area/sum(Area)) |> #why is it needed? Maybe can be removed
                        nest() |> 
                        mutate(models = map(data, ~lm(Area ~ `Analyte Concentration`, data = .x))) |> 
                        mutate(coef = map(models, coef)) |> 
                        mutate(Response_factor = map_dbl(models, ~ coef(.x)["`Analyte Concentration`"]))|> #get the slope
                        mutate(intercept = map(coef, pluck("(Intercept)"))) |> 
                        mutate(rsquared = map(models, summary)) |> 
                        mutate(rsquared = map(rsquared, pluck("r.squared"))) |>
                        select(-coef) |>  # remove coef variable since it has already been plucked
                        unnest(c(Response_factor, intercept, rsquared)) |>  #removing the list type for these variables
                        mutate(Response_factor = if_else(Response_factor < 0, 0, Response_factor)) |> # replace negative RF with 0
                        mutate(rsquared = ifelse(is.nan(rsquared), 0, rsquared)) |> 
                        #filter(rsquared >= removeRsquared()) |> 
                        mutate(Response_factor = if_else(rsquared < removeRsquared(), 0, Response_factor)) |>       
                        mutate(Chain_length = paste0("C", str_extract(Molecule, "(?<=C)[^H]+"))) |> 
                        #filter(Chain_length == "C10" | Chain_length == "C11" | Chain_length == "C12" | Chain_length == "C13") |> #this will be remove later or added as arg in fn
                        ungroup() |> 
                        group_by(!!input$standardAnnoColumn, Chain_length) |> #grouping by the selected Note
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
                
              
                
                CPs_samples <- Skyline_output() |>  # Call the reactive with ()
                        filter(`Sample Type` == "Unknown",
                               Molecule != "IS",
                               Molecule != "RS",
                               `Isotope Label Type` == "Quan") |> 
                        mutate(Chain_length = str_extract(Molecule, "(?<=C)[^H]+") |> as.numeric()) |>  # Extract number after "C"
                        #mutate(Group = case_when(
                         #       Chain_length >= 10 & Chain_length <= 13 ~ "S",  # Group S for 10-13
                         #       Chain_length >= 14 & Chain_length <= 17 ~ "M",  # Group M for 14-17
                          #      Chain_length >= 18 & Chain_length <= 30 ~ "L",  # Group L for 18-30
                           #     TRUE ~ "Unknown"  # Default case
                        #)) |> 
                        group_by(`Replicate Name`) |>  # Group by Replicate Name and Group
                        mutate(Relative_distribution = Area / sum(Area, na.rm = TRUE)) |> 
                        ungroup() |>  # Ungroup before dropping the Group column
                       #select(-Group) |>  # Explicitly remove Group column
                        select(`Replicate Name`, Molecule, Area, Relative_distribution) |>
                        mutate(across(Relative_distribution, ~replace(., is.nan(.), 0)))  # Replace NaN with zero
                
                
                CPs_standards_input <- CPs_standards |> 
                        select(Molecule, !!input$standardAnnoColumn, Response_factor) |> 
                        pivot_wider(names_from = !!input$standardAnnoColumn, values_from = "Response_factor")
                
                
                CPs_samples_input <- CPs_samples |> 
                        select(Molecule, `Replicate Name`, Relative_distribution) |> 
                        pivot_wider(names_from = "Replicate Name", values_from = "Relative_distribution")
                
                
                # This step ensures that all values are corresponding to the same molecule for std and sample        
                combined <- CPs_samples_input |> 
                        right_join(CPs_standards_input, by = "Molecule")
                
                ############################################################################### DECONVOLUTION #############################################################################
                
                # Ensure combined_matrix is correctly defined as a matrix prior to the deconvolution
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
                
                #Remove the replicate name to generate vectors:
                deconv_coef_df_matrix<- deconv_coef_df |> 
                        select(-`Replicate Name`)        
                
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
                df <- Skyline_output()  # Use the reactive data source
                req(df)  # Make sure df is not NULL
                
                # Prepare RECOVERY Data
                RECOVERY <- df |> 
                        filter(`Isotope Label Type` == "Quan") |> 
                        pivot_wider(id_cols = c(`Replicate Name`, `Sample Type`),
                                    names_from = Molecule, 
                                    values_from = Area) |>  # Spread IS and RS into columns
                        as.data.frame()  # Convert to data frame for consistency
                
                # Check if IS and RS columns exist before trying to mutate them
                if ("IS" %in% colnames(RECOVERY) & "RS" %in% colnames(RECOVERY)) {
                        RECOVERY <- RECOVERY |> 
                                mutate(across(c(IS, RS), ~replace_na(.x, 0)))  # Replace NAs with 0 for IS and RS
                } else {
                        stop("Required columns 'IS' or 'RS' are missing from the data.")
                }
                
                # Calculate AverageRatio for Quality Control
                dfR <- RECOVERY |> 
                        filter(`Sample Type` == "Quality Control") |> 
                        mutate(RatioStd = IS / RS) |> 
                        summarize(AverageRatio = mean(RatioStd, na.rm = TRUE))
                
                # Ensure dfR is valid
                if (nrow(dfR) == 0) {
                        stop("No Quality Control samples found.")
                }
                
                # RECOVERY Calculation
                RECOVERY <- RECOVERY |> 
                        filter(`Sample Type` %in% c("Unknown", "Blank")) |>  
                        mutate(RatioSample = IS / RS) |> 
                        mutate(Recovery = RatioSample / as.numeric(dfR$AverageRatio)) |>  
                        mutate(RecoveryPercentage = round(Recovery * 100, 0)) |>  # Round Recovery to 0 decimals
                        select(`Replicate Name`, `Sample Type`, RecoveryPercentage)  # Select only relevant columns
                
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
                        filter('Sample Type' == "Blank") |> 
                        summarize(LOD = sd(Concentration, na.rm = TRUE) * 3)  # Multiply by 3 for LOD
                
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

