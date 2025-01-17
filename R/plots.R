

plot_cal_SCCPs <- function(CPs_standards_SCCP) {
    # Prepare the data

    df <- CPs_standards_SCCP |>
        dplyr::filter(RF > 0) |>
        dplyr::mutate(fitted_values = purrr::map(models, purrr::pluck("fitted.values"))) |>
        tidyr::unnest(c(data, fitted_values)) |>
        dplyr::group_by(Batch_Name, Molecule) #need to group to get separate add_lines


    # Create the plot
    p <- plotly::plot_ly()

    # Add scatter plot markers
    p <- p |>
        plotly::add_markers(
            data = df,
            x = ~Analyte_Concentration,
            y = ~Area,
            color = ~PCA,
            #symbol = as.formula(paste0("~`", Batch_Name, "`")),
            symbol = ~Batch_Name,
            text = ~paste(
                "Homologue:", PCA,
                "<br>Area:", round(Area, 2),
                "<br>Analyte_Concentration:", round(Analyte_Concentration, 3),
                "<br>Standard:", Batch_Name,
                "<br>Rsquared:", round(rsquared, 3)
            ),
            hoverinfo = "text",
            visible = "legendonly",
            legendgroup = ~PCA
        )

    # Add linear model lines
    p <- p |>
        plotly::add_lines(
            data = df,
            x = ~Analyte_Concentration,
            y = ~fitted_values,
            color = ~PCA,
            line = list(dash = "solid"),
            hoverinfo = "none",
            visible = "legendonly",
            legendgroup = ~PCA,
            showlegend = FALSE  # Don't show duplicate legends for lines
        )

    # Configure layout
    p |>
        plotly::layout(
            title = list(
                text = "Calibration PCAs-C10-13",
                x = 0.5  # Center the title
            ),
            xaxis = list(
                title = "Analyte_Concentration",
                zeroline = TRUE,
                showgrid = TRUE
            ),
            yaxis = list(
                title = "Area",
                zeroline = TRUE,
                showgrid = TRUE
            )
        )
}



plot_cal_MCCPs <- function(CPs_standards_MCCP) {
    # Prepare the data

    df <- CPs_standards_MCCP |>
        dplyr::filter(RF > 0) |>
        dplyr::mutate(fitted_values = purrr::map(models, purrr::pluck("fitted.values"))) |>
        tidyr::unnest(c(data, fitted_values)) |>
        dplyr::group_by(Batch_Name, Molecule) #need to group to get separate add_lines


    # Create the plot
    p <- plotly::plot_ly()

    # Add scatter plot markers
    p <- p |>
        plotly::add_markers(
            data = df,
            x = ~Analyte_Concentration,
            y = ~Area,
            color = ~PCA,
            symbol = ~Batch_Name,
            text = ~paste(
                "Homologue:", PCA,
                "<br>Area:", round(Area, 2),
                "<br>Analyte_Concentration:", round(Analyte_Concentration, 3),
                "<br>Standard:", Batch_Name,
                "<br>Rsquared:", round(rsquared, 3)
            ),
            hoverinfo = "text",
            visible = "legendonly",
            legendgroup = ~PCA
        )

    # Add linear model lines
    p <- p |>
        plotly::add_lines(
            data = df,
            x = ~Analyte_Concentration,
            y = ~fitted_values,
            color = ~PCA,
            line = list(dash = "solid"),
            hoverinfo = "none",
            visible = "legendonly",
            legendgroup = ~PCA,
            showlegend = FALSE  # Don't show duplicate legends for lines
        )

    # Configure layout
    p |>
        plotly::layout(
            title = list(
                text = "Calibration PCAs-C14-17",
                x = 0.5  # Center the title
            ),
            xaxis = list(
                title = "Analyte_Concentration",
                zeroline = TRUE,
                showgrid = TRUE
            ),
            yaxis = list(
                title = "Area",
                zeroline = TRUE,
                showgrid = TRUE
            )
        )
}




plot_cal_LCCPs <- function(CPs_standards_LCCP) {
    # Prepare the data

    df <- CPs_standards_LCCP |>
        dplyr::filter(RF > 0) |>
        dplyr::mutate(fitted_values = purrr::map(models, purrr::pluck("fitted.values"))) |>
        tidyr::unnest(c(data, fitted_values)) |>
        dplyr::group_by(Batch_Name, Molecule) #need to group to get separate add_lines


    # Create the plot
    p <- plotly::plot_ly()

    # Add scatter plot markers
    p <- p |>
        plotly::add_markers(
            data = df,
            x = ~Analyte_Concentration,
            y = ~Area,
            color = ~PCA,
            symbol = ~Batch_Name,
            text = ~paste(
                "Homologue:", PCA,
                "<br>Area:", round(Area, 2),
                "<br>Analyte_Concentration:", round(Analyte_Concentration, 3),
                "<br>Standard:", Batch_Name,
                "<br>Rsquared:", round(rsquared, 3)
            ),
            hoverinfo = "text",
            visible = "legendonly",
            legendgroup = ~PCA
        )

    # Add linear model lines
    p <- p |>
        plotly::add_lines(
            data = df,
            x = ~Analyte_Concentration,
            y = ~fitted_values,
            color = ~PCA,
            line = list(dash = "solid"),
            hoverinfo = "none",
            visible = "legendonly",
            legendgroup = ~PCA,
            showlegend = FALSE  # Don't show duplicate legends for lines
        )

    # Configure layout
    p |>
        plotly::layout(
            title = list(
                text = "Calibration PCAs-C18-30",
                x = 0.5  # Center the title
            ),
            xaxis = list(
                title = "Analyte_Concentration",
                zeroline = TRUE,
                showgrid = TRUE
            ),
            yaxis = list(
                title = "Area",
                zeroline = TRUE,
                showgrid = TRUE
            )
        )
}
