# Various utilities and helper functions


####################################################################################
#-------------------------------- CPions functions --------------------------------#
####################################################################################

create_formula <- function(C, H, Cl, Br, S, O) {
    formula <- paste0(
        dplyr::case_when(C < 1 ~ paste0(""),
                  C == 1 ~ paste0("C"),
                  C > 1  ~ paste0("C", C)),
        dplyr::case_when(H < 1 ~ paste0(""),
                  H == 1 ~ paste0("H"),
                  H > 1  ~ paste0("H", H)),
        dplyr::case_when(Cl < 1 ~ paste0(""),
                  Cl == 1 ~ paste0("Cl"),
                  Cl > 1  ~ paste0("Cl", Cl)),
        dplyr::case_when(Br < 1 ~ paste0(""),
                  Br == 1 ~ paste0("Br"),
                  Br > 1  ~ paste0("Br", Br)),
        dplyr::case_when(S < 1 ~ paste0(""),
                  S == 1 ~ paste0("S"),
                  S > 1  ~ paste0("S", S)),
        dplyr::case_when(O < 1 ~ paste0(""),
                  O == 1 ~ paste0("O"),
                  O > 1  ~ paste0("O", O))
    )
    # Remove any leading or trailing spaces
    stringr::str_trim(formula)
}


create_elements <- function(data) {
    # String vector
    string_vector <- c("m/z", "abundance", "12C", "13C", "1H", "2H","35Cl", "37Cl", "79Br", "81Br", "16O", "17O", "18O", "32S", "33S", "34S", "36S")

    # Identify columns in the string vector that are not in the data frame
    new_columns <- base::setdiff(string_vector, names(data))

    # Add new columns to the data frame with default values: 0
    for (col in new_columns) {
        data[[col]] <- 0
    }
    return(data)
}

create_formula_isotope <- function(`12C`,`13C`, `1H`,`2H`, `35Cl`, `37Cl`, `79Br`, `81Br`, `16O`, `17O`, `18O`, `32S`, `33S`, `34S`, `36S`){
    formula_iso <- paste0(
        ifelse(`12C` > 0, paste0("[12C]", `12C`), ""),
        ifelse(`13C` > 0, paste0("[13C]", `13C`), ""),
        ifelse(`1H` > 0, paste0("[1H]", `1H`), ""),
        ifelse(`2H` > 0, paste0("[2H]", `2H`), ""),
        ifelse(`35Cl` > 0, paste0("[35Cl]", `35Cl`), ""),
        ifelse(`37Cl` > 0, paste0("[37Cl]", `37Cl`), ""),
        ifelse(`79Br` > 0, paste0("[79Br]", `79Br`), ""),
        ifelse(`81Br` > 0, paste0("[81Br]", `81Br`), ""),
        ifelse(`16O` > 0, paste0("[16O]", `16O`), ""),
        ifelse(`17O` > 0, paste0("[17O]", `17O`), ""),
        ifelse(`18O` > 0, paste0("[18O]", `18O`), ""),
        ifelse(`32S` > 0, paste0("[32S]", `32S`), ""),
        ifelse(`33S` > 0, paste0("[33S]", `33S`), ""),
        ifelse(`34S` > 0, paste0("[34S]", `34S`), ""),
        ifelse(`36S` > 0, paste0("[36S]", `36S`), "")
    )
    # Remove any leading or trailing spaces
    stringr::str_trim(formula_iso)
}


calculate_haloperc <- function(Molecule_Formula) {
    # Regular expression to extract atoms and their counts
    pattern <- "([A-Z][a-z]*)(\\d*)"
    mwtable <- data.frame(Atom = c("H", "C", "O", "S", "Cl", "Br"), MW = c(1.00794, 12.011, 15.9994, 32.066, 35.4527, 79.904))

    # Extract matches
    matches <- stringr::str_match_all(Molecule_Formula, pattern)[[1]]

    # Convert to a data frame for clarity
    result <- data.frame(
        Atom = matches[, 2],                     # Element symbols
        Count = as.numeric(matches[, 3])         # Element counts
    )

    # Replace missing counts (e.g., implicit "1") with 1
    result$Count[is.na(result$Count)] <- 1

    result <- result |>
        dplyr::left_join(mwtable, by = "Atom") |>
        dplyr::mutate(MW_atoms = MW*Count) |>
        dplyr::mutate(Halogen = case_when(Atom == "Cl" ~ TRUE,
                                   Atom == "Br" ~ TRUE,
                                   Atom == "F" ~ TRUE,
                                   Atom == "I" ~ TRUE,
                                   .default = FALSE))

    mw <- sum(result$MW_atoms)

    mw_halo <- result |>
        dplyr::filter(Halogen == TRUE) |>
        dplyr::summarise(sum(MW_atoms)) |>
        as.double()

    Molecule_Halo_perc <- round(mw_halo/mw*100, 0)
    return(Molecule_Halo_perc)
}


######### This function generates input for the Envipat function ###########

generateInput_Envipat_normal <- function(data = data, group = group, adduct_ions = adduct_ions, fragment_ions = fragment_ions) {


    data <- data |>
        dplyr::mutate(Halo_perc = dplyr::case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) |>
        dplyr::mutate(Adduct = adduct_ions) |>
        dplyr::mutate(Cl = dplyr::case_when(
            fragment_ions == "-Cl" ~ Cl-1,
            fragment_ions == "-HCl" ~ Cl-1,
            fragment_ions == "+Cl" ~ Cl+1,
            fragment_ions == "-Cl-HCl" ~ Cl-2,
            fragment_ions == "-Cl-2HCl" ~ Cl-3,
            fragment_ions == "-Cl-3HCl" ~ Cl-4,
            fragment_ions == "-Cl-4HCl" ~ Cl-5,
            fragment_ions == "-2Cl-HCl" ~ Cl-3,
            .default = Cl)) |>
        dplyr::mutate(H = dplyr::case_when(
            fragment_ions == "-H" ~ H-1,
            fragment_ions == "-HCl" ~ H-1,
            fragment_ions == "-Cl-HCl" ~ H-1,
            fragment_ions == "-Cl-2HCl" ~ H-2,
            fragment_ions == "-Cl-3HCl" ~ H-3,
            fragment_ions == "-Cl-4HCl" ~ H-4,
            fragment_ions == "-2Cl-HCl" ~ H-1,
            .default = H)) |>
        dplyr::mutate(Br = dplyr::case_when(
            fragment_ions == "+Br" ~ 1,
            TRUE ~0
        )) |>
        dplyr::mutate(Adduct_Formula = dplyr::case_when(
            fragment_ions != "+Br" ~ paste0("C", C, "H", H, "Cl", Cl),
            fragment_ions == "+Br" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br))) |>
        dplyr::select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, C, H, Cl)

    return(data)
}





# For bromochloroalkanes
generateInput_Envipat_BCA <- function(data = data, group = group, adduct_ions = adduct_ions, fragment_ions = fragment_ions) {


    data <- data |>
        dplyr::mutate(Halo_perc = round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0)) |>
        dplyr::mutate(Adduct = adduct_ions) |>
        dplyr::mutate(Cl = dplyr::case_when(
            fragment_ions == "-Cl" ~ Cl-1,
            fragment_ions == "-HCl" ~ Cl-1,
            fragment_ions == "+Cl" ~ Cl+1,
            fragment_ions == "-Cl-HCl" ~ Cl-2,
            fragment_ions == "-Cl-2HCl" ~ Cl-3,
            fragment_ions == "-Cl-3HCl" ~ Cl-4,
            fragment_ions == "-Cl-4HCl" ~ Cl-5,
            fragment_ions == "-2Cl-HCl" ~ Cl-3,
            .default = Cl)) |>
        dplyr::mutate(H = dplyr::case_when(
            fragment_ions == "-H" ~ H-1,
            fragment_ions == "-HCl" ~ H-1,
            fragment_ions == "-Cl-HCl" ~ H-1,
            fragment_ions == "-Cl-2HCl" ~ H-2,
            fragment_ions == "-Cl-3HCl" ~ H-3,
            fragment_ions == "-Cl-4HCl" ~ H-4,
            fragment_ions == "-2Cl-HCl" ~ H-1,
            .default = H)) |>
        dplyr::mutate(Adduct_Formula = paste0("C", C, "H", H, "Cl", Cl, "Br", Br)) |>
        dplyr::select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, C, H, Cl, Br)

    return(data)
}



generateInput_Envipat_advanced <- function(data = data, Compounds = Compounds, Adduct_Ion = Adduct_Ion,
                                           TP = TP, Charge = Charge) {


    data <- data |>
        dplyr::mutate(Adduct_Ion = Adduct_Ion) |>
        dplyr::mutate(Charge = Charge) |>
        dplyr::mutate(Adduct_Annotation = dplyr::case_when(
            TP == "None" ~ paste0("[", Compounds, Adduct_Ion, "]", Charge),
            .default = paste0("[", Compounds, TP, Adduct_Ion, "]", Charge))) |>
        dplyr::mutate(Adduct_Annotation = stringr::str_replace(Adduct_Annotation, "\\d$", "")) |>
        dplyr::mutate(Compound_Class = Compounds) |>
        dplyr::mutate(TP = TP) |>
        dplyr::mutate(Cl = dplyr::case_when(
            Adduct_Ion == "-Cl" ~ Cl-1,
            Adduct_Ion == "-HCl" ~ Cl-1,
            Adduct_Ion == "+Cl" ~ Cl+1,
            Adduct_Ion == "-Cl-HCl" ~ Cl-2,
            Adduct_Ion == "-Cl-2HCl" ~ Cl-3,
            Adduct_Ion == "-Cl-3HCl" ~ Cl-4,
            Adduct_Ion == "-Cl-4HCl" ~ Cl-5,
            Adduct_Ion == "-2Cl-HCl" ~ Cl-3,
            .default = Cl)) |>
        dplyr::mutate(H = dplyr::case_when(
            Adduct_Ion == "-H" ~ H-1,
            Adduct_Ion == "-HCl" ~ H-1,
            Adduct_Ion == "-Cl-HCl" ~ H-1,
            Adduct_Ion == "-Cl-2HCl" ~ H-2,
            Adduct_Ion == "-Cl-3HCl" ~ H-3,
            Adduct_Ion == "-Cl-4HCl" ~ H-4,
            Adduct_Ion == "-2Cl-HCl" ~ H-1,
            .default = H)) |>
        dplyr::mutate(Br = ifelse(Compound_Class == "BCA", Br, 0)) |>
        dplyr::mutate(Br = ifelse(Adduct_Ion == "+Br", Br+1, Br)) |>
        dplyr::mutate(O = dplyr::case_when(
            TP == "+OH" ~ 1,
            TP == "+2OH" ~ 2,
            TP == "+SO4" ~ 4,
            TP == "-Cl+OH" ~ 1,
            TP == "-2Cl+2OH" ~ 2,
            .default = 0
        )) |>
        dplyr::mutate(S = dplyr::case_when(
            TP == "+SO4" ~ 1,
            .default = 0)) |>
        dplyr::mutate(Adduct_Formula = create_formula(C, H, Cl, Br, S, O))|>
        dplyr::rowwise() |>
        dplyr::mutate(Molecule_Halo_perc = calculate_haloperc(Molecule_Formula)) |>
        dplyr::ungroup() |>
        dplyr::select(Molecule_Formula, Molecule_Halo_perc, Charge, Compound_Class, TP, Adduct_Ion, Adduct_Annotation, Adduct_Formula, C, H, Cl, Br, S, O)

    return(data)
}



getAdduct_normal <- function(adduct_ions, C, Cl, Clmax, threshold) {
    # Regex to extract strings
    ion_modes <- stringr::str_extract(adduct_ions, "(?<=\\]).{1}") # Using lookbehind assertion to extract ion mode
    fragment_ions <- stringr::str_extract(adduct_ions, "(?<=.{4}).+?(?=\\])") # extract after the 3rd character and before ]
    group <- stringr::str_extract(adduct_ions, "[^\\[].{2}") # Using positive lookbehind for [)

    if (group == "PCA") {
        data <- crossing(C, Cl) |> #set combinations of C and Cl
            dplyr::filter(C >= Cl) |> # filter so Cl dont exceed C atoms
            dplyr::filter(Cl <= Clmax) |> # limit chlorine atoms.
            dplyr::mutate(H = 2*C+2-Cl) |> # add H atoms
            dplyr::mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl)) |> #add chemical formula
            dplyr::select(Molecule_Formula, C, H, Cl) # move Formula to first column
    } else if (group == "PCO") {
        data <- crossing(C, Cl) |>
            dplyr::filter(C >= Cl) |>
            dplyr::filter(Cl <= Clmax) |>
            dplyr::mutate(H = 2*C-Cl) |>
            dplyr::mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl)) |>
            dplyr::select(Molecule_Formula, C, H, Cl)
    }  else {
        print("Input not correct, only PCA or PCO is allowed")
    }


    # adding ion modes to the data frame to be inserted to isopattern, only -1 or +1 allowed
    if (ion_modes == "-") {
        data <- data |>
            dplyr::mutate(Charge = as.integer(-1))
    }else if (ion_modes == "+") {
        data <- data |>
            dplyr::mutate(Charge = as.integer(1))
    }


    # generate input data for envipat based on fragment_ions
    data <- generateInput_Envipat_normal(data = data, group = group, adduct_ions = adduct_ions, fragment_ions = fragment_ions)

    # Remove formula without Cl after adduct formations
    data <- data |>
        dplyr::filter(Cl > 0)

    # Create empty list for all ion formulas
    CP_allions <- list()
    data_ls <- list()


    # function to get isotopic patterns for all PCAs.
    # data("isotopes") needs to be loaded in app.R
    getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes,
                                                    chemforms = x,
                                                    threshold = threshold,
                                                    emass = 0.00054857990924,
                                                    plotit = FALSE,
                                                    charge = Charge)}

    if (fragment_ions == "+Br") { #this is for Br adduct
        for (j in seq_along(data$Adduct_Formula)) {
            Adduct_Formula <- data$Adduct_Formula[j]
            Molecule_Formula <- data$Molecule_Formula[j]
            Charge <- data$Charge[j]
            Halo_perc <- data$Halo_perc[j]
            dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
            dat <- as.data.frame(dat[[1]])
            dat <- dat |>
                dplyr::mutate(abundance = round(abundance, 1)) |>
                dplyr::mutate(`m/z` = round(`m/z`, 6)) |>
                dplyr::mutate(Isotope_Formula = paste0("[12C]", `12C`, "[13C]", `13C`, "[1H]", `1H`, "[2H]", `2H`, "[35Cl]", `35Cl`, "[37Cl]", `37Cl`, "[79Br]", `79Br`, "[81Br]", `81Br`)) |>
                dplyr::mutate(Molecule_Formula = Molecule_Formula) |>
                dplyr::mutate(Halo_perc = Halo_perc) |>
                dplyr::mutate(Adduct_Formula =  Adduct_Formula) |>
                dplyr::mutate(Charge = Charge) |>
                dplyr::mutate(Isotopologue = dplyr::case_when(
                    `13C` + (`37Cl`+`81Br`)*2 == 0 ~ "",
                    `13C` + (`37Cl`+`81Br`)*2 == 1 ~ "+1",
                    `13C` + (`37Cl`+`81Br`)*2 == 2 ~ "+2",
                    `13C` + (`37Cl`+`81Br`)*2 == 3 ~ "+3",
                    `13C` + (`37Cl`+`81Br`)*2 == 4 ~ "+4",
                    `13C` + (`37Cl`+`81Br`)*2 == 5 ~ "+5",
                    `13C` + (`37Cl`+`81Br`)*2 == 6 ~ "+6",
                    `13C` + (`37Cl`+`81Br`)*2 == 7 ~ "+7",
                    `13C` + (`37Cl`+`81Br`)*2 == 8 ~ "+8",
                    `13C` + (`37Cl`+`81Br`)*2 == 9 ~ "+9",
                    `13C` + (`37Cl`+`81Br`)*2 == 10 ~ "+10",
                    `13C` + (`37Cl`+`81Br`)*2 == 11 ~ "+11",
                    `13C` + (`37Cl`+`81Br`)*2 == 12 ~ "+12",
                    `13C` + (`37Cl`+`81Br`)*2 == 13 ~ "+13",
                    `13C` + (`37Cl`+`81Br`)*2 == 14 ~ "+14",
                    `13C` + (`37Cl`+`81Br`)*2 == 15 ~ "+15",
                    `13C` + (`37Cl`+`81Br`)*2 == 16 ~ "+16",
                    `13C` + (`37Cl`+`81Br`)*2 == 17 ~ "+17",
                    `13C` + (`37Cl`+`81Br`)*2 == 18 ~ "+18",
                    `13C` + (`37Cl`+`81Br`)*2 == 19 ~ "+19",
                    `13C` + (`37Cl`+`81Br`)*2 == 20 ~ "+20")) |>
                dplyr::mutate(Adduct = paste0(adduct_ions, " ", Isotopologue)) |>
                dplyr::rename(Rel_ab = abundance) |>
                dplyr::select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`, `79Br`, `81Br`)
            data_ls[[j]] <- dat
        }
    }else { # for other adducts
        for (j in seq_along(data$Adduct_Formula)) {
            Adduct_Formula <- data$Adduct_Formula[j]
            Molecule_Formula <- data$Molecule_Formula[j]
            Charge <- data$Charge[j]
            Halo_perc <- data$Halo_perc[j]
            dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
            dat <- as.data.frame(dat[[1]])
            dat <- dat |>
                dplyr::mutate(abundance = round(abundance, 1)) |>
                dplyr::mutate(`m/z` = round(`m/z`, 6)) |>
                dplyr::mutate(Isotope_Formula = paste0("[12C]", `12C`, "[13C]", `13C`, "[1H]", `1H`, "[2H]", `2H`, "[35Cl]", `35Cl`, "[37Cl]", `37Cl`)) |>
                dplyr::mutate(Molecule_Formula = Molecule_Formula) |>
                dplyr::mutate(Halo_perc = Halo_perc) |>
                dplyr::mutate(Adduct_Formula =  Adduct_Formula) |>
                dplyr::mutate(Charge = Charge) |>
                dplyr::mutate(Isotopologue = dplyr::case_when(
                    `13C` + (`37Cl`)*2 == 0 ~ "",
                    `13C` + (`37Cl`)*2 == 1 ~ "+1",
                    `13C` + (`37Cl`)*2 == 2 ~ "+2",
                    `13C` + (`37Cl`)*2 == 3 ~ "+3",
                    `13C` + (`37Cl`)*2 == 4 ~ "+4",
                    `13C` + (`37Cl`)*2 == 5 ~ "+5",
                    `13C` + (`37Cl`)*2 == 6 ~ "+6",
                    `13C` + (`37Cl`)*2 == 7 ~ "+7",
                    `13C` + (`37Cl`)*2 == 8 ~ "+8",
                    `13C` + (`37Cl`)*2 == 9 ~ "+9",
                    `13C` + (`37Cl`)*2 == 10 ~ "+10",
                    `13C` + (`37Cl`)*2 == 11 ~ "+11",
                    `13C` + (`37Cl`)*2 == 12 ~ "+12",
                    `13C` + (`37Cl`)*2 == 13 ~ "+13",
                    `13C` + (`37Cl`)*2 == 14 ~ "+14",
                    `13C` + (`37Cl`)*2 == 15 ~ "+15",
                    `13C` + (`37Cl`)*2 == 16 ~ "+16",
                    `13C` + (`37Cl`)*2 == 17 ~ "+17",
                    `13C` + (`37Cl`)*2 == 18 ~ "+18",
                    `13C` + (`37Cl`)*2 == 19 ~ "+19",
                    `13C` + (`37Cl`)*2 == 20 ~ "+20")) |>
                dplyr::mutate(Adduct = paste0(adduct_ions, " ", Isotopologue)) |>
                dplyr::rename(Rel_ab = abundance) |>
                dplyr::select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`)
            data_ls[[j]] <- dat
        }
    }

    # combine all elements in list list to get dataframe
    data_ls <- do.call(rbind, data_ls)


    # combine both all adduct ions
    CP_allions <- rbind(CP_allions, data_ls)
    return(CP_allions)

}


getAdduct_BCA <- function(adduct_ions, C, Cl, Br, Clmax, Brmax, threshold) {

    # Regex to extract strings
    ion_modes <- stringr::str_extract(adduct_ions, "(?<=\\]).{1}") # Using lookbehind assertion to extract ion mode
    fragment_ions <- stringr::str_extract(adduct_ions, "(?<=.{4}).+?(?=\\])") # extract after the 3rd character and before ]
    group <- stringr::str_extract(adduct_ions, "[^\\[].{2}") # Using positive lookbehind for [)

    if (group == "BCA") {
        data <- crossing(C, Cl, Br) |> #get combinations of C, Cl, Br
            dplyr::filter(C >= Cl) |> # filter so Cl dont exceed C atoms
            dplyr::filter(Cl <= Clmax) |> # limit chlorine atoms.
            dplyr::filter(Br <= Brmax) |>
            dplyr::filter(Br + Cl <= C) |>
            dplyr::mutate(H = 2*C+2-Cl-Br) |> # add H atoms
            dplyr::mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl, "Br", Br)) |> #add chemical formula
            dplyr::select(Molecule_Formula, C, H, Cl, Br) # move Formula to first column
    }  else {
        print("Input not correct, only BCA is allowed")
    }

    # check chem_forms
    # if (any(check_chemform(isotopes = isotopes, chemforms = data$Formula)$warning == TRUE)) {print("Warning: incorrect formula")} else {"All correct"}

    # adding ion modes to the data frame to be inserted to isopattern, only -1 or +1 allowed
    if (ion_modes == "-") {
        data <- data |>
            dplyr::mutate(Charge = as.integer(-1))
    }else if (ion_modes == "+") {
        data <- data |>
            dplyr::mutate(Charge = as.integer(1))
    }


    ####### generate input data for envipat based on fragment_ions
    data <- generateInput_Envipat_BCA(data = data, group = group, adduct_ions = adduct_ions, fragment_ions = fragment_ions)



    # Remove formula without Cl after adduct formations
    data <- data |>
        dplyr::filter(Cl > 0)

    # Create empty list for all ion formulas
    CP_allions <- list()
    data_ls <- list()


    # function to get isotopic patterns for all PCAs. Threshold based on the app, neutral form. data("isotopes") needs to be loaded first
    getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes,
                                                    chemforms = x,
                                                    threshold = threshold,
                                                    emass = 0.00054857990924,
                                                    plotit = FALSE,
                                                    charge = Charge)}


    for (j in seq_along(data$Adduct_Formula)) {
        Adduct_Formula <- data$Adduct_Formula[j]
        Molecule_Formula <- data$Molecule_Formula[j]
        Charge <- data$Charge[j]
        Halo_perc <- data$Halo_perc[j]
        dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
        dat <- as.data.frame(dat[[1]])
        dat <- dat |>
            dplyr::mutate(abundance = round(abundance, 1)) |>
            dplyr::mutate(`m/z` = round(`m/z`, 6)) |>
            dplyr::mutate(Isotope_Formula = paste0("[12C]", `12C`, "[13C]", `13C`, "[1H]", `1H`, "[2H]", `2H`, "[35Cl]", `35Cl`, "[37Cl]", `37Cl`, "[79Br]", `79Br`, "[81Br]", `81Br`)) |>
            dplyr::mutate(Molecule_Formula = Molecule_Formula) |>
            dplyr::mutate(Halo_perc = Halo_perc) |>
            dplyr::mutate(Adduct_Formula = Adduct_Formula) |>
            dplyr::mutate(Charge = Charge) |>
            dplyr::mutate(Isotopologue = dplyr::case_when(
                `13C` + (`37Cl`+`81Br`)*2 == 0 ~ "",
                `13C` + (`37Cl`+`81Br`)*2 == 1 ~ "+1",
                `13C` + (`37Cl`+`81Br`)*2 == 2 ~ "+2",
                `13C` + (`37Cl`+`81Br`)*2 == 3 ~ "+3",
                `13C` + (`37Cl`+`81Br`)*2 == 4 ~ "+4",
                `13C` + (`37Cl`+`81Br`)*2 == 5 ~ "+5",
                `13C` + (`37Cl`+`81Br`)*2 == 6 ~ "+6",
                `13C` + (`37Cl`+`81Br`)*2 == 7 ~ "+7",
                `13C` + (`37Cl`+`81Br`)*2 == 8 ~ "+8",
                `13C` + (`37Cl`+`81Br`)*2 == 9 ~ "+9",
                `13C` + (`37Cl`+`81Br`)*2 == 10 ~ "+10",
                `13C` + (`37Cl`+`81Br`)*2 == 11 ~ "+11",
                `13C` + (`37Cl`+`81Br`)*2 == 12 ~ "+12",
                `13C` + (`37Cl`+`81Br`)*2 == 13 ~ "+13",
                `13C` + (`37Cl`+`81Br`)*2 == 14 ~ "+14",
                `13C` + (`37Cl`+`81Br`)*2 == 15 ~ "+15",
                `13C` + (`37Cl`+`81Br`)*2 == 16 ~ "+16",
                `13C` + (`37Cl`+`81Br`)*2 == 17 ~ "+17",
                `13C` + (`37Cl`+`81Br`)*2 == 18 ~ "+18",
                `13C` + (`37Cl`+`81Br`)*2 == 19 ~ "+19",
                `13C` + (`37Cl`+`81Br`)*2 == 20 ~ "+20")) |>
            dplyr::mutate(Adduct = paste0(adduct_ions, " ", Isotopologue)) |>
            dplyr::rename(Rel_ab = abundance) |>
            dplyr::select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`, `79Br`, `81Br`)
        data_ls[[j]] <- dat
    }


    # combine all elements in list list to get dataframe
    data_ls <- do.call(rbind, data_ls)


    # combine both all adduct ions
    CP_allions <- rbind(CP_allions, data_ls)
    return(CP_allions)

}



getAdduct_advanced <- function(Compounds, Adduct_Ion, TP, Charge, C, Cl, Clmax, Br, Brmax, threshold) {

    # Regex to extract strings
    if (Compounds == "PCA") {
        data <- crossing(C, Cl) |> #set combinations of C and Cl
            dplyr::filter(C >= Cl) |> # filter so Cl dont exceed C atoms
            dplyr::filter(Cl <= Clmax) |> # limit chlorine atoms.
            dplyr::mutate(H = dplyr::case_when(# add H atoms
                TP == "None" ~ 2*C+2-Cl,
                TP == "+OH" ~ 2*C+2-Cl,
                TP == "+2OH" ~ 2*C+2-Cl,
                TP == "-Cl+OH" ~ 2*C+2-Cl+1,
                TP == "-2Cl+2OH" ~ 2*C+2-Cl+2,
                TP == "+SO4" ~ 2*C+2-Cl-1))  |>
            dplyr::mutate(Cl = dplyr::case_when(
                TP == "-Cl+OH" ~ Cl-1,
                TP == "-2Cl+2OH" ~ Cl-2,
                .default = Cl)) |>
            dplyr::mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl)) |>
            dplyr::mutate(Molecule_Formula = case_when( #DOUBLE CHECK THE FORMULA IS CORRECT!!!!!
                TP == "None" ~ paste0("C", C, "H", H, "Cl", Cl),
                TP == "+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                TP == "+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                TP == "-Cl+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                TP == "-2Cl+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                TP == "+SO4" ~ paste0("C", C, "H", H, "Cl", Cl, "SO4")))

    } else if (Compounds == "PCO") {
        data <- crossing(C, Cl) |>
            dplyr::filter(C >= Cl) |>
            dplyr::filter(Cl <= Clmax) |>
            dplyr::mutate(H = dplyr::case_when(# add H atoms.
                TP == "None" ~ 2*C-Cl, #PCO general formula 2*C-Cl
                TP == "+OH" ~ 2*C-Cl,
                TP == "+2OH" ~ 2*C-Cl,
                TP == "-Cl+OH" ~ 2*C-Cl+1,
                TP == "-2Cl+2OH" ~ 2*C-Cl+2,
                TP == "+SO4" ~ 2*C-Cl-1))  |>
            dplyr::mutate(Cl = dplyr::case_when(
                TP == "-Cl+OH" ~ Cl-1,
                TP == "-2Cl+2OH" ~ Cl-2,
                .default = Cl)) |>
            dplyr::mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl)) |>
            dplyr::mutate(Molecule_Formula = dplyr::case_when( #DOUBLE CHECK THE FORMULA IS CORRECT!!!!!
                TP == "None" ~ paste0("C", C, "H", H, "Cl", Cl),
                TP == "+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                TP == "+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                TP == "-Cl+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                TP == "-2Cl+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                TP == "+SO4" ~ paste0("C", C, "H", H, "Cl", Cl, "SO4")))

    } else if (Compounds == "BCA") {
        data <- tidyr::crossing(C, Cl, Br) |>  #get combinations of C, Cl, Br
            dplyr::filter(C >= Cl) |>  # filter so Cl dont exceed C atoms
            dplyr::filter(Cl <= Clmax) |>  # limit chlorine atoms.
            dplyr::filter(Br <= Brmax) |>
            dplyr::filter(Br + Cl <= C) |>
            dplyr::mutate(H = dplyr::case_when(# add H atoms.
                TP == "None" ~ 2*C+2-Cl-Br, #BCA general formula
                TP == "+OH" ~ 2*C+2-Cl-Br,
                TP == "+2OH" ~ 2*C+2-Cl-Br,
                TP == "-Cl+OH" ~ 2*C+2-Cl-Br+1,
                TP == "-2Cl+2OH" ~ 2*C+2-Cl-Br+2,
                TP == "+SO4" ~ 2*C+2-Cl-Br-1))  |>
            dplyr::mutate(Molecule_Formula = dplyr::case_when( #DOUBLE CHECK THE FORMULA IS CORRECT!!!!!
                TP == "None" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br),
                TP == "+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O"),
                TP == "+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O2"),
                TP == "-Cl+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O"),
                TP == "-2Cl+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O2"),
                TP == "+SO4" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "SO4")))


    }



    # adding ion modes to the data frame to be inserted to isopattern, only -1 or +1 allowed
    if (Charge == "-") {
        data <- data |>
            dplyr::mutate(Charge = as.integer(-1))
    }else if (Charge == "+") {
        data <- data |>
            dplyr::mutate(Charge = as.integer(1))
    }


    #generate input data for envipat based on adduct ions
    data <- generateInput_Envipat_advanced(data = data, Compounds = Compounds, Adduct_Ion = Adduct_Ion, TP = TP, Charge = Charge)

    # Remove formula without Cl after adduct formations
    data <- data |>
        dplyr::filter(Cl > 0)

    # Create empty list for all ion formulas
    CP_allions <- list()
    data_ls <- list()


    # function to get isotopic patterns for all PCA/PCO/BCA.
    # data("isotopes") needs to be loaded in app.R
    getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes,
                                                    chemforms = x,
                                                    threshold = threshold,
                                                    emass = 0.00054857990924,
                                                    plotit = FALSE,
                                                    charge = Charge)}

    for (j in seq_along(data$Molecule_Formula)) {
        Adduct_Formula <- data$Adduct_Formula[j]
        Molecule_Formula <- data$Molecule_Formula[j]
        Compound_Class <- data$Compound_Class[j]
        TP <- data$TP[j]
        Charge <- data$Charge[j]
        Molecule_Halo_perc <- data$Molecule_Halo_perc[j]
        Adduct_Annotation <- data$Adduct_Annotation[j]
        dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
        dat <- as.data.frame(dat[[1]])

        dat <- dat |>
            dplyr::mutate(abundance = round(abundance, 1)) |>
            dplyr::mutate(`m/z` = round(`m/z`, 6))

        dat <- create_elements(dat) |>
            dplyr::mutate(Isotope_Formula = create_formula_isotope(`12C`,`13C`, `1H`, `2H`, `35Cl`, `37Cl`, `79Br`, `81Br`,
                                                            `16O`, `17O`, `18O`, `32S`, `33S`, `34S`, `36S`)) |>
            dplyr::mutate(Molecule_Formula = Molecule_Formula) |>
            dplyr::mutate(Molecule_Halo_perc = Molecule_Halo_perc) |>
            dplyr::mutate(Compound_Class = Compound_Class) |>
            dplyr::mutate(TP = TP) |>
            dplyr::mutate(Adduct_Annotation =  Adduct_Annotation) |>
            dplyr::mutate(Adduct_Formula =  Adduct_Formula) |>
            dplyr::mutate(Charge = Charge) |>
            dplyr::mutate(Isotopologue = case_when(
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 0 ~ "",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 1 ~ "+1",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 2 ~ "+2",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 3 ~ "+3",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 4 ~ "+4",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 5 ~ "+5",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 6 ~ "+6",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 7 ~ "+7",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 8 ~ "+8",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 9 ~ "+9",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 10 ~ "+10",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 11 ~ "+11",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 12 ~ "+12",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 13 ~ "+13",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 14 ~ "+14",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 15 ~ "+15",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 16 ~ "+16",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 17 ~ "+17",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 18 ~ "+18",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 19 ~ "+19",
                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 20 ~ "+20")) |>
            dplyr::mutate(Adduct_Isotopologue = paste0(Adduct_Ion, " ", Isotopologue)) |>
            dplyr::rename(Rel_ab = abundance) |>
            dplyr::select(Molecule_Formula, Molecule_Halo_perc, Compound_Class, TP, Charge, Adduct_Annotation, Adduct_Isotopologue, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`, everything())
        data_ls[[j]] <- dat
    }


    # combine all elements in list list to get dataframe
    data_ls <- do.call(rbind, data_ls)


    # combine both all adduct ions
    CP_allions <- rbind(CP_allions, data_ls)
    return(CP_allions)

}


#########################################################################################################
#------------------------------------------ CPquant functions ------------------------------------------#
#########################################################################################################

### Function to perform deconvolution on a single data frame ###
perform_deconvolution <- function(df, combined_standard) {
    df_matrix <- as.matrix(df)


    print(paste("df_matrix dimensions:", dim(df_matrix)))
    print(paste("combined_standard dimensions:", dim(combined_standard)))

    if (nrow(combined_standard) != nrow(df_matrix)) {
        stop("Dimensions of combined_standard and df are incompatible.")
    }

    # Reshape df_matrix if it has only one column or extract the first column if it has multiple
    if (ncol(df_matrix) == 1) {
        df_vector <- as.vector(df_matrix)
    } else {
        df_vector <- as.vector(df_matrix["Relative_distribution"])  # Extract the first column for nnls
    }

    # Check for NA/NaN/Inf values in df_vector and combined_standard
    if (any(is.na(df_vector)) || any(is.nan(df_vector)) || any(is.infinite(df_vector))) {
        stop("df_vector contains NA/NaN/Inf values.")
    }

    if (any(is.na(combined_standard)) || any(is.nan(combined_standard)) || any(is.infinite(combined_standard))) {
        stop("combined_standard contains NA/NaN/Inf values.")
    }

    # Perform nnls
    deconv <- nnls::nnls(combined_standard, df_vector)

    # Extract deconvolution results
    deconv_coef <- deconv$x

    # # Normalize the coefficients so they sum to 100%
    # if (sum(deconv_coef) > 0) {
    #     deconv_coef <- deconv_coef / sum(deconv_coef) * 100
    # }

    # Normalize the coefficients so they sum to 1
    if (sum(deconv_coef) > 0) {
     deconv_coef <- deconv_coef / sum(deconv_coef)
    }

    # Calculate deconvolved values using matrix multiplication
    deconv_resolved <- combined_standard %*% deconv_coef



    # Calculate the Goodnes of fit
    # Calculate the total sum of squares (SST)
    sst <- sum((df_vector - mean(df_vector))^2)
    # Calculate the residual sum of squares (SSR)
    ssr <- sum((df_vector - deconv_resolved/sum(deconv_resolved))^2)

    # Calculate R-squared goodnes of fit
    r_squared <- 1 - (ssr / sst)
    # Calculate Mean Squared Error (MSE)
    mse <- mean((df_vector - deconv_resolved/sum(deconv_resolved))^2)

    # Calculate Root Mean Squared Error (RMSE)
    rmse <- sqrt(mse)

    # Chiq-square test: ensure that values are positive for chi-square test
    # if (any(deconv_resolved < 0) || any(df_vector < 0)) {
    #     warning("Non-positive values found, skipping chi-square test")
    #     chisq_result <- NULL
    # } else {
    #     #adding a small constant to avoid 0 values
    #     observed_corr <- df_vector + 1E-9
    #     predicted_corr <- deconv_resolved + 1E-9
    #     chisq_result <- chisq.test(x= observed_corr, p = predicted_corr/sum(predicted_corr), rescale.p = TRUE)
    # }

    combined_standard_names <- colnames(combined_standard)

    names(deconv_coef) <- combined_standard_names


    return(list(
        deconv_coef = deconv_coef,
        deconv_resolved = deconv_resolved,
        r_squared = r_squared
        #chisq_result = chisq_result
    ))
}

