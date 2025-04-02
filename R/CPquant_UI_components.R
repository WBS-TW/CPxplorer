## CPquant UI Components ##
###########################

defineVariablesUI <- function(Skyline_output){
###START: Define UI components

    # Create the UI components
    shiny::fluidRow(
        shiny::column(
            6,
            shiny::selectInput(
                inputId = "removeSamples", #select if some samples will be removed from quantification
                label = 'Remove samples from quantification?',
                choices = unique(Skyline_output$Replicate_Name),
                selected = NULL,
                multiple = TRUE
            )
        ),

        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
        shiny::column(
            6,
            shiny::sliderInput(
                inputId = "removeRsquared", #keep only Molecule above this rsquared, zero means keep everything
                label = 'Keep the the calibration curves above this rsquared (0 means keep everything)',
                min = 0,
                max = 1,
                value = 0.80,
                step = 0.05
            )
        )
    )
}
### END FUNCTION


defineRecoveryUI <- function(Skyline_output){
    ###START: Define UI components

    # Create the UI components
    shiny::fluidRow(
        shiny::column(
            6,
            shiny::selectInput(
                inputId = "chooseRS", #select which will be the RS
                label = 'Choose RS',
                choices = unique(Skyline_output$Molecule[Skyline_output$Molecule_List == "RS"]),
                selected = NULL,
                multiple = FALSE
            )
        ))
}
### END FUNCTION

