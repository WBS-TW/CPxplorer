
#############################################################################
plot_skyline_output <- function(Skyline_output){

    Skyline_output |>
        dplyr::filter(Isotope_Label_Type == "Quan") |>
        dplyr::mutate(OrderedMolecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)]))) |>  # Create a composite ordering factor
        plotly::plot_ly(
            x = ~ OrderedMolecule,
            y = ~ Area,
            color = ~ Sample_Type,
            type = "box",
            text = ~paste(
                "Homologue: ", PCA,
                "<br>Sample: ", Replicate_Name,
                "<br>Area:", round(Area, 2)
            ),
            hoverinfo = "text"
        ) |>
        plotly::layout(xaxis = list(title = 'Molecule'),
               yaxis = list(title = 'Area'))
}

#############################################################################

plot_calibration_curves <- function(CPs_standards, quantUnit) {
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

        plot_list[[i]] <- plotly::plot_ly() |>
            plotly::add_trace(
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
                    "<br>R2:", round(rsquared, 3)
                ),
                hoverinfo = 'text'
            ) |>
            plotly::add_trace(
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
        shareX = FALSE,
        #shareY = TRUE,
        shareY = FALSE,
        margin = 0.1
    ) |>
        plotly::layout(
            height = 400 * subplot_rows,
            showlegend = TRUE,
            annotations = annotations,
            margin = list(t = 50, b = 50, l = 50, r = 50),
            legend = list(
                groupclick = "togglegroup",  # Changed from "toggleitem" to "togglegroup"
                tracegroupgap = 10,
                itemsizing = "constant"
            ),
            xaxis = list(title = paste0("Analyte Concentration/Amount (", quantUnit, ")"))  # Use the reactive value for x-axis label
        )

    return(final_plot)
}

#############################################################################

plot_quanqualratio <- function(Skyline_output_filt) {

Skyline_output_filt |>
    dplyr::group_by(Replicate_Name, Molecule) |>
    dplyr::mutate(Quan_Area = ifelse(Isotope_Label_Type == "Quan", Area, NA)) |>
    tidyr::fill(Quan_Area, .direction = "downup") |>
    dplyr::mutate(QuanMZ = ifelse(Isotope_Label_Type == "Quan", Chromatogram_Precursor_MZ, NA)) |>
    tidyr::fill(QuanMZ, .direction = "downup") |>
    dplyr::mutate(QuanQualRatio = ifelse(Isotope_Label_Type == "Qual", Quan_Area/Area, 1)) |>
    tidyr::replace_na(list(QuanQualRatio = 0)) |>
    dplyr::mutate(QuanQualMZ = paste0(QuanMZ,"/",Chromatogram_Precursor_MZ)) |>
    dplyr::ungroup() |>
    dplyr::select(Replicate_Name, Sample_Type, Molecule_List, Molecule, QuanQualMZ, QuanQualRatio) |>
    plotly::plot_ly(x = ~Replicate_Name, y = ~QuanQualRatio, type = 'violin', color = ~Sample_Type,
            text = ~paste("Sample: ", Replicate_Name,
                          "<br>Molecule List: ", Molecule_List,
                          "<br>Molecule: ", Molecule,
                          "<br>Quan/Qual MZ: ", QuanQualMZ,
                          "<br>Ratio: ", round(QuanQualRatio, 2)),
            hoverinfo = "text") |>
    plotly::layout(title = 'Quan-to-Qual Ratio',
           xaxis = list(title = 'Replicate Name'),
           yaxis = list(title = 'Quan-to-Qual Ratio'))
}

##############################################################################

plot_meas_vs_theor_ratio <- function(Skyline_output_filt) {

    Skyline_output_filt |>
        dplyr::group_by(Replicate_Name, Molecule) |>
        dplyr::mutate(Quan_Area = ifelse(Isotope_Label_Type == "Quan", Area, NA)) |>
        tidyr::fill(Quan_Area, .direction = "downup") |>
        dplyr::mutate(QuanMZ = ifelse(Isotope_Label_Type == "Quan", Chromatogram_Precursor_MZ, NA)) |>
        tidyr::fill(QuanMZ, .direction = "downup") |>
        dplyr::mutate(QuanQualRatio = ifelse(Isotope_Label_Type == "Qual", Quan_Area/Area, 1)) |>
        tidyr::replace_na(list(QuanQualRatio = 0)) |>
        dplyr::mutate(QuanQualMZ = paste0(QuanMZ,"/",Chromatogram_Precursor_MZ)) |>

        dplyr::mutate(Quan_Rel_Ab = ifelse(Isotope_Label_Type == "Quan", Rel_Ab, NA)) |>
        tidyr::fill(Quan_Rel_Ab, .direction = "downup") |>
        dplyr::mutate(QuanQual_Rel_Ab_Ratio = ifelse(Isotope_Label_Type == "Qual", Quan_Rel_Ab/Rel_Ab, 1)) |>
        tidyr::replace_na(list(QuanQual_Rel_Ab_Ratio = 0)) |>
        dplyr::ungroup() |>
        dplyr::mutate(MeasVSTheo = QuanQualRatio/QuanQual_Rel_Ab_Ratio) |>
        dplyr::mutate(Is_Outlier = MeasVSTheo > 3 | MeasVSTheo < 0.3) |>
        dplyr::mutate(Is_Outlier = factor(Is_Outlier, levels = c(FALSE, TRUE), labels = c("Within Limit", "Outlier"))) |>
        dplyr::select(Replicate_Name, Sample_Type, Molecule_List, Molecule, QuanQualMZ, QuanQualRatio, QuanQual_Rel_Ab_Ratio, MeasVSTheo, Is_Outlier) |>
        plotly::plot_ly(x = ~Replicate_Name, y = ~MeasVSTheo,
                        type = 'scatter', mode = 'markers',
                        color = ~Is_Outlier,
                        colors = c('blue', 'red'),
                        text = ~paste("Replicate:", Replicate_Name,
                                      "<br>Homologue Group: ", Molecule,
                                      "<br>Measured against Theoretical Ratio:", round(MeasVSTheo, 1)),
                        marker = list(size = 10)) |>
        layout(title = "Measured/Theoretical ratio >3 or <0.3 are marked in red (ratio of 1 means perfect match",
               xaxis = list(title = "Sample Name"),
               yaxis = list(title = "MeasVSTheo"))

}


##############################################################################

plot_sample_contribution <- function(deconvolution) {

    # How much contribution of each sample to the final deconvoluted homologue group pattern
    plot_data <- deconvolution |>
        tidyr::unnest(deconv_coef) |>
        tidyr::unnest_longer(c(deconv_coef, Batch_Name)) |>
        dplyr::select(Replicate_Name, Batch_Name, deconv_coef)

    # Create the plotly stacked bar plot
    plotly::plot_ly(plot_data,
            x = ~Replicate_Name,
            y = ~deconv_coef,
            type = "bar",
            color = ~Batch_Name,
            colors = "Spectral") |>
        plotly::layout(
            title = list(
                text = "Contributions from standards to deconvoluted homologue pattern",
                x = 0.5,  # Center the title
                y = 0.95  # Position slightly down from top
            ),
            barmode = "stack",
            xaxis = list(title = "Replicate Name"),
            yaxis = list(title = "Relative Contribution",
                         tickformat = ".2%"),
            showlegend = TRUE,
            legend = list(title = list(text = "Batch Name"))
        )
}

#############################################################################

plot_homologue_group_pattern_comparison <- function(Sample_distribution, input_selectedSamples){

    # Filter data for selected samples and reshape data
    selected_samples <- Sample_distribution |>
        dplyr::filter(Replicate_Name %in% input_selectedSamples) |>
        dplyr::mutate(Molecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)])))

    # Get unique homologue groups for consistent coloring
    homologue_groups <- unique(selected_samples$C_homologue)

    # Create a list of plots, one for each Replicate_Name
    plot_list <- selected_samples |>
        split(selected_samples$Replicate_Name) |>
        map(function(df) {
            # Create base plot
            p <- plotly::plot_ly()

            # Add bars for each homologue group
            for(hg in homologue_groups) {
                df_filtered <- df[df$C_homologue == hg,]
                p <- p |>
                    plotly::add_trace(
                        data = df_filtered,
                        x = ~Molecule,
                        y = ~Relative_Area,
                        name = hg,
                        legendgroup = hg,
                        showlegend = (df_filtered$Replicate_Name[1] == input_selectedSamples[1]),
                        type = 'bar',
                        opacity = 1
                    )
            }

            # Add the black line for resolved_distribution
            p <- p |>
                plotly::add_trace(
                    data = df,
                    x = ~Molecule,
                    y = ~resolved_distribution,
                    name = "Deconvoluted Distribution",
                    legendgroup = "DeconvDistr",
                    showlegend = (df$Replicate_Name[1] == input_selectedSamples[1]),
                    type = 'scatter',
                    mode = 'lines+markers',
                    line = list(color = 'black'),
                    marker = list(color = 'black', size = 6),
                    opacity = 0.7
                )

            # Add layout
            p <- p |>
                plotly::layout(
                    xaxis = list(
                        title = "Homologue",
                        tickangle = 45
                    ),
                    yaxis = list(title = "Value"),
                    barmode = 'group',
                    annotations = list(
                        x = 0.5,
                        y = 1.1,
                        text = unique(df$Replicate_Name),
                        xref = 'paper',
                        yref = 'paper',
                        showarrow = FALSE
                    )
                )

            return(p)
        })

    # Combine the plots using subplot
    plotly::subplot(plot_list,
                    nrows = ceiling(length(plot_list)/2),
                    shareX = TRUE,
                    shareY = TRUE) |>
        plotly::layout(
            #title = "Sample Comparison",
            showlegend = TRUE,
            hovermode = 'closest',
            hoverlabel = list(bgcolor = "white"),
            barmode = 'group'
        ) |>
        plotly::config(displayModeBar = TRUE) |>
        htmlwidgets::onRender("
                function(el) {
                    var plotDiv = document.getElementById(el.id);
                    plotDiv.on('plotly_legendclick', function(data) {
                        Plotly.restyle(plotDiv, {
                            visible: data.data[data.curveNumber].visible === 'legendonly' ? true : 'legendonly'
                        }, data.fullData.map((trace, i) => i).filter(i =>
                            data.fullData[i].legendgroup === data.fullData[data.curveNumber].legendgroup
                        ));
                        return false;
                    });
                }
            ")
}

