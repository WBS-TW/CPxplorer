

plot_cal_SCCPs <- function(CPs_standards_S, standardAnnoColumn) {
    # Prepare the data

    df <- CPs_standards_S |>
        dplyr::filter(Response_factor > 0) |>
        #dplyr::filter(`Molecule List` %in% c("PCA-C10", "PCA-C11", "PCA-C12", "PCA-C13")) |>
        dplyr::mutate(fitted_values = purrr::map(models, purrr::pluck("fitted.values"))) |>
        tidyr::unnest(c(data, fitted_values)) |>
        dplyr::group_by(!!dplyr::sym(standardAnnoColumn), `Molecule`) #need to group to get separate add_lines

    # Ensure standardAnnoColumn exists in the data
    if (!standardAnnoColumn %in% colnames(df)) {
        stop("standardAnnoColumn not found in data")
    }

    # Create the plot
    p <- plotly::plot_ly()

    # Add scatter plot markers
    p <- p |>
        plotly::add_markers(
            data = df,
            x = ~`Analyte Concentration`,
            y = ~Area,
            color = ~PCA,
            symbol = as.formula(paste0("~`", standardAnnoColumn, "`")),
            text = ~paste(
                "Homologue:", PCA,
                "<br>Area:", round(Area, 2),
                "<br>Analyte Concentration:", round(`Analyte Concentration`, 3),
                "<br>Standard:", get(standardAnnoColumn),
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
            x = ~`Analyte Concentration`,
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
                title = "Analyte Concentration",
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



plot_cal_MCCPs <- function(CPs_standards_M, standardAnnoColumn) {
    # Prepare the data

    df <- CPs_standards_M |>
        dplyr::filter(Response_factor > 0) |>
        #dplyr::filter(`Molecule List` %in% c("PCA-C14", "PCA-C15", "PCA-C16", "PCA-C17")) |>
        dplyr::mutate(fitted_values = purrr::map(models, purrr::pluck("fitted.values"))) |>
        tidyr::unnest(c(data, fitted_values)) |>
        dplyr::group_by(!!dplyr::sym(standardAnnoColumn), `Molecule`) #need to group to get separate add_lines

    # Ensure standardAnnoColumn exists in the data
    if (!standardAnnoColumn %in% colnames(df)) {
        stop("standardAnnoColumn not found in data")
    }

    # Create the plot
    p <- plotly::plot_ly()

    # Add scatter plot markers
    p <- p |>
        plotly::add_markers(
            data = df,
            x = ~`Analyte Concentration`,
            y = ~Area,
            color = ~PCA,
            symbol = as.formula(paste0("~`", standardAnnoColumn, "`")),
            text = ~paste(
                "Homologue:", PCA,
                "<br>Area:", round(Area, 2),
                "<br>Analyte Concentration:", round(`Analyte Concentration`, 3),
                "<br>Standard:", get(standardAnnoColumn),
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
            x = ~`Analyte Concentration`,
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
                title = "Analyte Concentration",
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




plot_cal_LCCPs <- function(CPs_standards_L, standardAnnoColumn) {
    # Prepare the data

    df <- CPs_standards_L |>
        dplyr::filter(Response_factor > 0) |>
        # dplyr::filter(`Molecule List` %in% c("PCA-C18", "PCA-C19", "PCA-C20", "PCA-C21", "PCA-C22", "PCA-C23",
        #                                      "PCA-C24", "PCA-C25", "PCA-C26", "PCA-C27", "PCA-C28", "PCA-C29", "PCA-C30")) |>
        dplyr::mutate(fitted_values = purrr::map(models, purrr::pluck("fitted.values"))) |>
        tidyr::unnest(c(data, fitted_values)) |>
        dplyr::group_by(!!dplyr::sym(standardAnnoColumn), `Molecule`) #need to group to get separate add_lines

    # Ensure standardAnnoColumn exists in the data
    if (!standardAnnoColumn %in% colnames(df)) {
        stop("standardAnnoColumn not found in data")
    }

    # Create the plot
    p <- plotly::plot_ly()

    # Add scatter plot markers
    p <- p |>
        plotly::add_markers(
            data = df,
            x = ~`Analyte Concentration`,
            y = ~Area,
            color = ~PCA,
            symbol = as.formula(paste0("~`", standardAnnoColumn, "`")),
            text = ~paste(
                "Homologue:", PCA,
                "<br>Area:", round(Area, 2),
                "<br>Analyte Concentration:", round(`Analyte Concentration`, 3),
                "<br>Standard:", get(standardAnnoColumn),
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
            x = ~`Analyte Concentration`,
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
                title = "Analyte Concentration",
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
