
plot_calibration_curves <- function(CPs_standards) {
    # Unnest the data first
    CPs_standards_unnested <- CPs_standards |>
        dplyr::filter(RF > 0) |>
        tidyr::unnest(data) |>
        dplyr::mutate(Molecule = factor(Molecule,
                                        levels = unique(Molecule[order(C_number, Cl_number)])))

    # Get unique Quantification Groups
    groups <- unique(CPs_standards_unnested$Quantification_Group)

    # Create a list to store individual plots
    plot_list <- list()

    # Create individual plots for each group
    for(i in seq_along(groups)) {
        group_data <- CPs_standards_unnested |>
            dplyr::filter(Quantification_Group == groups[i])

        plot_list[[i]] <- plot_ly() |>
            add_trace(
                data = group_data,
                x = ~Analyte_Concentration,
                y = ~Area,
                color = ~PCA,
                type = 'scatter',
                mode = 'markers',
                legendgroup = ~Molecule_List,
                legendgrouptitle = list(text = ~Molecule_List),
                showlegend = FALSE,
                name = ~paste(PCA, "(points)"),
                text = ~paste(
                    "Molecule:", Molecule,
                    "<br>Area:", round(Area, 2),
                    "<br>Concentration:", round(Analyte_Concentration, 2),
                    "<br>R²:", round(rsquared, 3)
                ),
                hoverinfo = 'text'
            ) |>
            add_trace(
                data = group_data,
                x = ~Analyte_Concentration,
                y = ~RF * Analyte_Concentration + intercept,
                color = ~PCA,
                type = 'scatter',
                mode = 'lines',
                legendgroup = ~Molecule_List,
                name = ~PCA,
                hoverinfo = 'none'
            )
    }

    # Calculate layout
    subplot_cols <- min(length(groups), 2)  # Maximum 2 columns
    subplot_rows <- ceiling(length(groups) / subplot_cols)

    # Create annotations for titles
    annotations <- list()
    for(i in seq_along(groups)) {
        row <- ceiling(i/subplot_cols)
        col <- if(i %% subplot_cols == 0) subplot_cols else i %% subplot_cols

        annotations[[i]] <- list(
            text = groups[i],
            font = list(size = 14),
            xref = "paper",
            yref = "paper",
            x = (col - 0.5)/subplot_cols,
            y = 1 - (row - 1)/subplot_rows,
            xanchor = "center",
            yanchor = "bottom",
            showarrow = FALSE
        )
    }

    # Combine plots using subplot
    final_plot <- plotly::subplot(
        plot_list,
        nrows = subplot_rows,
        shareX = TRUE,
        shareY = TRUE,
        margin = 0.1
    ) |>
        layout(
            height = 400 * subplot_rows,
            showlegend = TRUE,
            annotations = annotations,
            margin = list(t = 50, b = 50, l = 50, r = 50),
            legend = list(
                groupclick = "togglegroup",  # Changed from "toggleitem" to "togglegroup"
                tracegroupgap = 10,
                itemsizing = "constant"
            )
        )

    return(final_plot)
}




