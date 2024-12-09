
load_plots <- function(Skyline_output) {


    output$CalibrationSCCPs <- plotly::renderPlotly({
        data <- Skyline_output() |>
            dplyr::filter(`Isotope Label Type` == "Quan") |>
            dplyr::filter(stringr::str_detect(Note, "S-")) |>
            dplyr::filter(`Molecule List` %in% c("PCA-C10", "PCA-C11", "PCA-C12", "PCA-C13")) |>
            dplyr::group_by(Note, `Molecule`)

        # Fit linear models for each combination of Note and Molecule List
        models <- data |>
            dplyr::group_modify(~ {
                lm_model <- stats::lm(Area ~ `Analyte Concentration`, data = .x)
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
                lm_model <- stats::lm(Area ~ `Analyte Concentration`, data = .x)
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

}
