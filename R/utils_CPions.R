# Various utilities and helper functions for CPions



####################################################################################
#-------------------------------- CPions functions --------------------------------#
####################################################################################


#############################################################################
#############################################################################


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


#############################################################################
#############################################################################


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


#############################################################################
#############################################################################


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


#############################################################################
#############################################################################


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



#############################################################################
#############################################################################
######### This function generates input for the Envipat function ###########


generateInput_Envipat_normal <- function(data = data, group = group, adduct_ions = adduct_ions, fragment_ions = fragment_ions) {

    data <- data |>
        dplyr::mutate(Halo_perc = dplyr::case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) |>
        dplyr::mutate(Compound_Class = dplyr::case_when(group == "PCA" ~ "PCA",
                                                        group == "PCO" ~"PCO",
                                                        group == "BCA" ~ "BCA")) |>
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
        dplyr::select(Molecule_Formula, Compound_Class, Halo_perc, Charge, Adduct, Adduct_Formula, C, H, Cl)

    return(data)
}




#############################################################################
#############################################################################


generateInput_Envipat_BCA <- function(data = data, group = group, adduct_ions = adduct_ions, fragment_ions = fragment_ions) {


    data <- data |>
        dplyr::mutate(Halo_perc = round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0)) |>
        dplyr::mutate(Compound_Class = dplyr::case_when(group == "PCA" ~ "PCA",
                                                        group == "PCO" ~"PCO",
                                                        group == "BCA" ~ "BCA")) |>
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
        dplyr::select(Molecule_Formula, Compound_Class, Halo_perc, Charge, Adduct, Adduct_Formula, C, H, Cl, Br)

    return(data)
}


#############################################################################
#############################################################################


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
            TP == "-H+OH" ~ 1,
            TP == "-2H+2OH" ~ 2,
            TP == "-H+SO4H" ~ 4,
            TP == "-Cl+OH" ~ 1,
            TP == "-2Cl+2OH" ~ 2,
            TP == "-2H+O" ~ 1,
            .default = 0
        )) |>
        dplyr::mutate(S = dplyr::case_when(
            TP == "-H+SO4H" ~ 1,
            .default = 0)) |>
        dplyr::mutate(Adduct_Formula = create_formula(C, H, Cl, Br, S, O))|>
        dplyr::rowwise() |>
        dplyr::mutate(Molecule_Halo_perc = calculate_haloperc(Molecule_Formula)) |>
        dplyr::ungroup() |>
        dplyr::select(Molecule_Formula, Molecule_Halo_perc, Charge, Compound_Class, TP, Adduct_Ion, Adduct_Annotation, Adduct_Formula, C, H, Cl, Br, S, O)

    return(data)
}


#############################################################################
#############################################################################


getAdduct_normal <- function(adduct_ions, C, Cl, Clmax, threshold) {
    # Regex to extract strings
    ion_modes <- stringr::str_extract(adduct_ions, "(?<=\\]).{1}") # Using lookbehind assertion to extract ion mode
    fragment_ions <- stringr::str_extract(adduct_ions, "(?<=.{4}).+?(?=\\])") # extract after the 3rd character and before ]
    group <- stringr::str_extract(adduct_ions, "(?<=\\[)[A-Za-z]+(?=[+-])") # Using positive lookbehind precedes a [ ; matches on or more letters ; positive lookahead of either + or -

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
            Compound_Class <- data$Compound_Class[j]
            Charge <- data$Charge[j]
            Halo_perc <- data$Halo_perc[j]
            dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
            dat <- as.data.frame(dat[[1]])
            dat <- dat |>
                dplyr::mutate(Compound_Class = Compound_Class) |>
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
                dplyr::select(Molecule_Formula, Compound_Class, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`, `79Br`, `81Br`)
            data_ls[[j]] <- dat
        }
    }else { # for other adducts
        for (j in seq_along(data$Adduct_Formula)) {
            Adduct_Formula <- data$Adduct_Formula[j]
            Molecule_Formula <- data$Molecule_Formula[j]
            Compound_Class <- data$Compound_Class[j]
            Charge <- data$Charge[j]
            Halo_perc <- data$Halo_perc[j]
            dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
            dat <- as.data.frame(dat[[1]])
            dat <- dat |>
                dplyr::mutate(Compound_Class = Compound_Class) |>
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
                dplyr::select(Molecule_Formula, Compound_Class, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`)
            data_ls[[j]] <- dat
        }
    }

    # combine all elements in list list to get dataframe
    data_ls <- do.call(rbind, data_ls)


    # combine both all adduct ions
    CP_allions <- rbind(CP_allions, data_ls)
    return(CP_allions)

}


#############################################################################
#############################################################################


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
        Compound_Class <- data$Compound_Class[j]
        Adduct_Formula <- data$Adduct_Formula[j]
        Molecule_Formula <- data$Molecule_Formula[j]
        Charge <- data$Charge[j]
        Halo_perc <- data$Halo_perc[j]
        dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
        dat <- as.data.frame(dat[[1]])
        dat <- dat |>
            dplyr::mutate(Compound_Class = Compound_Class) |>
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
            dplyr::select(Molecule_Formula, Compound_Class, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`, `79Br`, `81Br`)
        data_ls[[j]] <- dat
    }


    # combine all elements in list list to get dataframe
    data_ls <- do.call(rbind, data_ls)


    # combine both all adduct ions
    CP_allions <- rbind(CP_allions, data_ls)
    return(CP_allions)

}

#############################################################################
#############################################################################


getAdduct_advanced <- function(Compounds, Adduct_Ion, TP, Charge, C, Cl, Clmax, Br, Brmax, threshold) {

    # Regex to extract strings
    if (Compounds == "PCA") {
        data <- crossing(C, Cl) |> #set combinations of C and Cl
            dplyr::filter(C >= Cl) |> # filter so Cl dont exceed C atoms
            dplyr::filter(Cl <= Clmax) |> # limit chlorine atoms.
            dplyr::mutate(H = dplyr::case_when(# add H atoms
                TP == "None" ~ 2*C+2-Cl, #PCA general formula 2*C+2-Cl
                TP == "-H+OH" ~ 2*C+2-Cl,
                TP == "-2H+2OH" ~ 2*C+2-Cl,
                TP == "-Cl+OH" ~ 2*C+2-Cl+1,
                TP == "-2Cl+2OH" ~ 2*C+2-Cl+2,
                TP == "-2H+O" ~ 2*C-Cl,
                TP == "-H+SO4H" ~ 2*C+2-Cl))  |>
            dplyr::mutate(Cl = dplyr::case_when(
                TP == "-Cl+OH" ~ Cl-1,
                TP == "-2Cl+2OH" ~ Cl-2,
                .default = Cl)) |>
            dplyr::mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl)) |>
            dplyr::mutate(Molecule_Formula = case_when(
                TP == "None" ~ paste0("C", C, "H", H, "Cl", Cl),
                TP == "-H+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                TP == "-2H+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                TP == "-Cl+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                TP == "-2Cl+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                TP == "-2H+O" ~ paste0("C", C, "H", H, "Cl", Cl,"O"),
                TP == "-H+SO4H" ~ paste0("C", C, "H", H, "Cl", Cl, "SO4")))

    } else if (Compounds == "PCO") {
        data <- crossing(C, Cl) |>
            dplyr::filter(C >= Cl) |>
            dplyr::filter(Cl <= Clmax) |>
            dplyr::mutate(H = dplyr::case_when(# add H atoms.
                TP == "None" ~ 2*C-Cl, #PCO general formula 2*C-Cl
                TP == "-H+OH" ~ 2*C-Cl,
                TP == "-2H+2OH" ~ 2*C-Cl,
                TP == "-Cl+OH" ~ 2*C-Cl+1,
                TP == "-2Cl+2OH" ~ 2*C-Cl+2,
                TP == "-2H+O" ~ 2*C-Cl-2,
                TP == "-H+SO4H" ~ 2*C-Cl))  |>
            dplyr::mutate(Cl = dplyr::case_when(
                TP == "-Cl+OH" ~ Cl-1,
                TP == "-2Cl+2OH" ~ Cl-2,
                .default = Cl)) |>
            dplyr::mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl)) |>
            dplyr::mutate(Molecule_Formula = dplyr::case_when( #DOUBLE CHECK THE FORMULA IS CORRECT!!!!!
                TP == "None" ~ paste0("C", C, "H", H, "Cl", Cl),
                TP == "-H+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                TP == "-2H+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                TP == "-Cl+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                TP == "-2Cl+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                TP == "-2H+O" ~ paste0("C", C, "H", H, "Cl", Cl,"O"),
                TP == "-H+SO4H" ~ paste0("C", C, "H", H, "Cl", Cl, "SO4")))

    } else if (Compounds == "BCA") {
        data <- tidyr::crossing(C, Cl, Br) |>  #get combinations of C, Cl, Br
            dplyr::filter(C >= Cl) |>  # filter so Cl dont exceed C atoms
            dplyr::filter(Cl <= Clmax) |>  # limit chlorine atoms.
            dplyr::filter(Br <= Brmax) |>
            dplyr::filter(Br + Cl <= C) |>
            dplyr::mutate(H = dplyr::case_when(# add H atoms.
                TP == "None" ~ 2*C+2-Cl-Br, #BCA general formula
                TP == "-H+OH" ~ 2*C+2-Cl-Br,
                TP == "-2H+2OH" ~ 2*C+2-Cl-Br,
                TP == "-Cl+OH" ~ 2*C+2-Cl-Br+1,
                TP == "-2Cl+2OH" ~ 2*C+2-Cl-Br+2,
                TP == "-2H+O" ~ 2*C-Cl-Br,
                TP == "-H+SO4H" ~ 2*C+2-Cl-Br))  |>
            dplyr::mutate(Molecule_Formula = dplyr::case_when( #DOUBLE CHECK THE FORMULA IS CORRECT!!!!!
                TP == "None" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br),
                TP == "-H+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O"),
                TP == "-2H+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O2"),
                TP == "-Cl+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O"),
                TP == "-2Cl+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O2"),
                TP == "-2H+O" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O"),
                TP == "-H+SO4H" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "SO4")))


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

########################################################################
########################################################################

# from envipat
isotopes <- structure(list(element = c("H", "H", "He", "He", "Li", "Li",
                           "Be", "B", "B", "C", "C", "N", "N", "O", "O", "O", "F", "Ne",
                           "Ne", "Ne", "Na", "Mg", "Mg", "Mg", "Al", "Si", "Si", "Si", "P",
                           "S", "S", "S", "S", "S", "Cl", "Cl", "Cl", "Ar", "Ar", "Ar",
                           "K", "K", "K", "Ca", "Ca", "Ca", "Ca", "Ca", "Ca", "Ca", "Ca",
                           "Ca", "Sc", "Ti", "Ti", "Ti", "Ti", "Ti", "V", "V", "Cr", "Cr",
                           "Cr", "Cr", "Mn", "Fe", "Fe", "Fe", "Fe", "Fe", "Co", "Ni", "Ni",
                           "Ni", "Ni", "Ni", "Cu", "Cu", "Zn", "Zn", "Zn", "Zn", "Zn", "Ga",
                           "Ga", "Ge", "Ge", "Ge", "Ge", "Ge", "As", "Se", "Se", "Se", "Se",
                           "Se", "Se", "Br", "Br", "Br", "Kr", "Kr", "Kr", "Kr", "Kr", "Kr",
                           "Rb", "Rb", "Sr", "Sr", "Sr", "Sr", "Y", "Zr", "Zr", "Zr", "Zr",
                           "Zr", "Nb", "Mo", "Mo", "Mo", "Mo", "Mo", "Mo", "Mo", "Ru", "Ru",
                           "Ru", "Ru", "Ru", "Ru", "Ru", "Rh", "Pd", "Pd", "Pd", "Pd", "Pd",
                           "Pd", "Ag", "Ag", "Cd", "Cd", "Cd", "Cd", "Cd", "Cd", "Cd", "Cd",
                           "In", "In", "Sn", "Sn", "Sn", "Sn", "Sn", "Sn", "Sn", "Sn", "Sn",
                           "Sn", "Sb", "Sb", "Te", "Te", "Te", "Te", "Te", "Te", "Te", "Te",
                           "I", "Xe", "Xe", "Xe", "Xe", "Xe", "Xe", "Xe", "Xe", "Xe", "Cs",
                           "Ba", "Ba", "Ba", "Ba", "Ba", "Ba", "Ba", "La", "La", "Ce", "Ce",
                           "Ce", "Ce", "Pr", "Nd", "Nd", "Nd", "Nd", "Nd", "Nd", "Nd", "Sm",
                           "Sm", "Sm", "Sm", "Sm", "Sm", "Sm", "Eu", "Eu", "Gd", "Gd", "Gd",
                           "Gd", "Gd", "Gd", "Gd", "Tb", "Dy", "Dy", "Dy", "Dy", "Dy", "Dy",
                           "Dy", "Ho", "Er", "Er", "Er", "Er", "Er", "Er", "Tm", "Yb", "Yb",
                           "Yb", "Yb", "Yb", "Yb", "Yb", "Lu", "Lu", "Hf", "Hf", "Hf", "Hf",
                           "Hf", "Hf", "Ta", "Ta", "W", "W", "W", "W", "W", "Re", "Re",
                           "Os", "Os", "Os", "Os", "Os", "Os", "Os", "Ir", "Ir", "Pt", "Pt",
                           "Pt", "Pt", "Pt", "Pt", "Au", "Hg", "Hg", "Hg", "Hg", "Hg", "Hg",
                           "Hg", "Tl", "Tl", "Pb", "Pb", "Pb", "Pb", "Bi", "Th", "Pa", "U",
                           "U", "U", "[15]N", "[13]C", "D", "[37]Cl", "[18]O", "[16]O",
                           "[12]C", "[2]H", "[33]S", "[34]S", "[35]S", "[36]S", "[35]Cl"
), isotope = c("1H", "2H", "3He", "4He", "6Li", "7Li", "9Be",
               "10B", "11B", "12C", "13C", "14N", "15N", "16O", "17O", "18O",
               "19F", "20Ne", "21Ne", "22Ne", "23Na", "24Mg", "25Mg", "26Mg",
               "27Al", "28Si", "29Si", "30Si", "31P", "32S", "33S", "34S", "35S",
               "36S", "35Cl", "36Cl", "37Cl", "36Ar", "38Ar", "40Ar", "39K",
               "40K", "41K", "40Ca", "41Ca", "42Ca", "43Ca", "44Ca", "45Ca",
               "46Ca", "47Ca", "48Ca", "45Sc", "46Ti", "47Ti", "48Ti", "49Ti",
               "50Ti", "50V", "51V", "50Cr", "52Cr", "53Cr", "54Cr", "55Mn",
               "54Fe", "55Fe", "56Fe", "57Fe", "58Fe", "59Co", "58Ni", "60Ni",
               "61Ni", "62Ni", "64Ni", "63Cu", "65Cu", "64Zn", "66Zn", "67Zn",
               "68Zn", "70Zn", "69Ga", "71Ga", "70Ge", "72Ge", "73Ge", "74Ge",
               "76Ge", "75As", "74Se", "76Se", "77Se", "78Se", "80Se", "82Se",
               "79Br", "80Br", "81Br", "78Kr", "80Kr", "82Kr", "83Kr", "84Kr",
               "86Kr", "85Rb", "87Rb", "84Sr", "86Sr", "87Sr", "88Sr", "89Y",
               "90Zr", "91Zr", "92Zr", "94Zr", "96Zr", "93Nb", "92Mo", "94Mo",
               "95Mo", "96Mo", "97Mo", "98Mo", "100Mo", "96Ru", "98Ru", "99Ru",
               "100Ru", "101Ru", "102Ru", "104Ru", "103Rh", "102Pd", "104Pd",
               "105Pd", "106Pd", "108Pd", "110Pd", "107Ag", "109Ag", "106Cd",
               "108Cd", "110Cd", "111Cd", "112Cd", "113Cd", "114Cd", "116Cd",
               "113In", "115In", "112Sn", "114Sn", "115Sn", "116Sn", "117Sn",
               "118Sn", "119Sn", "120Sn", "122Sn", "124Sn", "121Sb", "123Sb",
               "120Te", "122Te", "123Te", "124Te", "125Te", "126Te", "128Te",
               "130Te", "127I", "124Xe", "126Xe", "128Xe", "129Xe", "130Xe",
               "131Xe", "132Xe", "134Xe", "136Xe", "133Cs", "130Ba", "132Ba",
               "134Ba", "135Ba", "136Ba", "137Ba", "138Ba", "138La", "139La",
               "136Ce", "138Ce", "140Ce", "142Ce", "141Pr", "142Nd", "143Nd",
               "144Nd", "145Nd", "146Nd", "148Nd", "150Nd", "144Sm", "147Sm",
               "148Sm", "149Sm", "150Sm", "152Sm", "154Sm", "151Eu", "153Eu",
               "152Gd", "154Gd", "155Gd", "156Gd", "157Gd", "158Gd", "160Gd",
               "159Tb", "156Dy", "158Dy", "160Dy", "161Dy", "162Dy", "163Dy",
               "164Dy", "165Ho", "162Er", "164Er", "166Er", "167Er", "168Er",
               "170Er", "169Tm", "168Yb", "170Yb", "171Yb", "172Yb", "173Yb",
               "174Yb", "176Yb", "175Lu", "176Lu", "174Hf", "176Hf", "177Hf",
               "178Hf", "179Hf", "180Hf", "180Ta", "181Ta", "180W", "182W",
               "183W", "184W", "186W", "185Re", "187Re", "184Os", "186Os", "187Os",
               "188Os", "189Os", "190Os", "192Os", "191Ir", "193Ir", "190Pt",
               "192Pt", "194Pt", "195Pt", "196Pt", "198Pt", "197Au", "196Hg",
               "198Hg", "199Hg", "200Hg", "201Hg", "202Hg", "204Hg", "203Tl",
               "205Tl", "204Pb", "206Pb", "207Pb", "208Pb", "209Bi", "232Th",
               "231Pa", "234U", "235U", "238U", "15N", "13C", "2H", "37Cl",
               "18O", "16O", "12C", "2H", "33S", "34S", "35S", "36S", "35Cl"
), mass = c(1.007825032, 2.014101778, 3.016029319, 4.002603254,
            6.0151223, 7.0160041, 9.0121822, 10.0129371, 11.0093055, 12,
            13.00335484, 14.00307401, 15.00010897, 15.99491462, 16.9991315,
            17.9991604, 18.9984032, 19.99244018, 20.99384668, 21.99138511,
            22.98976966, 23.98504187, 24.985837, 25.982593, 26.98153863,
            27.97692649, 28.97649468, 29.97377018, 30.97376149, 31.97207073,
            32.97145854, 33.96786687, 35, 35.96708088, 34.96885271, 36, 36.9659026,
            35.96754511, 37.9627324, 39.96238312, 38.9637069, 39.96399867,
            40.96182597, 39.9625912, 41, 41.9586183, 42.9587668, 43.9554811,
            45, 45.9536927, 47, 47.952533, 44.9559119, 45.9526316, 46.9517631,
            47.9479463, 48.94787, 49.9447912, 49.9471585, 50.9439595, 49.9460442,
            51.9405075, 52.9406494, 53.9388804, 54.9380451, 53.9396147, 55,
            55.9349418, 56.9353983, 57.9332801, 58.933195, 57.9353429, 59.9307864,
            60.931056, 61.9283451, 63.927966, 62.9295975, 64.9277895, 63.9291422,
            65.9260334, 66.9271273, 67.9248442, 69.9253193, 68.9255736, 70.9247013,
            69.9242474, 71.9220758, 72.9234589, 73.9211778, 75.9214026, 74.9215965,
            73.9224764, 75.9192136, 76.919914, 77.9173091, 79.9165213, 81.9166994,
            78.9183379, 80, 80.916291, 77.9203648, 79.916379, 81.9134836,
            82.914136, 83.911507, 85.91061073, 84.91178974, 86.90918053,
            83.913425, 85.9092602, 86.9088771, 87.9056121, 88.9058483, 89.9047044,
            90.9056458, 91.9050408, 93.9063152, 95.9082734, 92.9063781, 91.906811,
            93.9050883, 94.9058421, 95.9046795, 96.9060215, 97.9054082, 99.90747,
            95.907598, 97.905287, 98.9059393, 99.9042195, 100.9055821, 101.9043493,
            103.905433, 102.905504, 101.905609, 103.904036, 104.905085, 105.903486,
            107.903892, 109.905153, 106.905097, 108.904752, 105.906459, 107.904184,
            109.9030021, 110.9041781, 111.9027578, 112.9044017, 113.9033585,
            115.904756, 112.904058, 114.903878, 111.904818, 113.902779, 114.903342,
            115.901741, 116.902952, 117.901603, 118.903308, 119.9021947,
            121.903439, 123.9052739, 120.9038157, 122.904214, 119.90402,
            121.9030439, 122.90427, 123.9028179, 124.9044307, 125.9033117,
            127.9044631, 129.9062244, 126.904473, 123.905893, 125.904274,
            127.9035313, 128.9047794, 129.903508, 130.9050824, 131.9041535,
            133.9053945, 135.907219, 132.9054519, 129.9063208, 131.9050613,
            133.9045084, 134.9056886, 135.9045759, 136.9058274, 137.9052472,
            137.907112, 138.9063533, 135.907172, 137.905991, 139.9054387,
            141.909244, 140.9076528, 141.9077233, 142.9098143, 143.9100873,
            144.9125736, 145.9131169, 147.916893, 149.920891, 143.911999,
            146.9148979, 147.9148227, 148.9171847, 149.9172755, 151.9197324,
            153.9222093, 150.9198502, 152.9212303, 151.919791, 153.9208656,
            154.922622, 155.9221227, 156.9239601, 157.9241039, 159.9270541,
            158.9253468, 155.924283, 157.924409, 159.9251975, 160.9269334,
            161.9267984, 162.9287312, 163.9291748, 164.9303221, 161.928778,
            163.9292, 165.9302931, 166.9320482, 167.9323702, 169.9354643,
            168.9342133, 167.933897, 169.9347618, 170.9363258, 171.9363815,
            172.9382108, 173.9388621, 175.9425717, 174.9407718, 175.9426863,
            173.940046, 175.9414086, 176.9432207, 177.9436988, 178.9458161,
            179.94655, 179.9474648, 180.9479958, 179.946704, 181.9482042,
            182.950223, 183.9509312, 185.9543641, 184.952955, 186.9557531,
            183.9524891, 185.9538382, 186.9557505, 187.9558382, 188.9581475,
            189.958447, 191.9614807, 190.960594, 192.9629264, 189.959932,
            191.961038, 193.9626803, 194.9647911, 195.9649515, 197.967893,
            196.9665687, 195.965833, 197.966769, 198.9682799, 199.968326,
            200.9703023, 201.970643, 203.9734939, 202.9723442, 204.9744275,
            203.9730436, 205.9744653, 206.9758969, 207.9766521, 208.9803987,
            232.0380553, 231.035884, 234.0409521, 235.0439299, 238.0507882,
            15.000109, 13.003355, 2.01410178, 36.9659026, 17.9991604, 15.99491462,
            12, 2.014101778, 32.97145854, 33.96786687, 35, 35.96708088, 34.96885271
), abundance = c(0.999885, 0.000115, 1.34e-06, 0.99999866, 0.0759,
                 0.9241, 1, 0.199, 0.801, 0.9893, 0.0107, 0.99636, 0.00364, 0.99757,
                 0.00038, 0.00205, 1, 0.9048, 0.0027, 0.0925, 1, 0.7899, 0.1,
                 0.1101, 1, 0.92223, 0.04685, 0.03092, 1, 0.9499, 0.0075, 0.0425,
                 0, 1e-04, 0.7576, 0, 0.2424, 0.003365, 0.000632, 0.996003, 0.932581,
                 0.000117, 0.067302, 0.96941, 0, 0.00647, 0.00135, 0.02086, 0,
                 4e-05, 0, 0.00187, 1, 0.0825, 0.0744, 0.7372, 0.0541, 0.0518,
                 0.0025, 0.9975, 0.04345, 0.83789, 0.09501, 0.02365, 1, 0.05845,
                 0, 0.91754, 0.02119, 0.00282, 1, 0.680769, 0.262231, 0.011399,
                 0.036345, 0.009256, 0.6915, 0.3085, 0.48268, 0.27975, 0.04102,
                 0.19024, 0.00631, 0.60108, 0.39892, 0.2038, 0.2731, 0.0776, 0.3672,
                 0.0783, 1, 0.0089, 0.0937, 0.0763, 0.2377, 0.4961, 0.0873, 0.5069,
                 0, 0.4931, 0.00355, 0.02286, 0.11593, 0.115, 0.56987, 0.17279,
                 0.7217, 0.2783, 0.0056, 0.0986, 0.07, 0.8258, 1, 0.5145, 0.1122,
                 0.1715, 0.1738, 0.028, 1, 0.1477, 0.0923, 0.159, 0.1668, 0.0956,
                 0.2419, 0.0967, 0.0554, 0.0187, 0.1276, 0.126, 0.1706, 0.3155,
                 0.1862, 1, 0.0102, 0.1114, 0.2233, 0.2733, 0.2646, 0.1172, 0.51839,
                 0.48161, 0.0125, 0.0089, 0.1249, 0.128, 0.2413, 0.1222, 0.2873,
                 0.0749, 0.0429, 0.9571, 0.0097, 0.0066, 0.0034, 0.1454, 0.0768,
                 0.2422, 0.0859, 0.3258, 0.0463, 0.0579, 0.5721, 0.4279, 9e-04,
                 0.0255, 0.0089, 0.0474, 0.0707, 0.1884, 0.3174, 0.3408, 1, 0.000952,
                 0.00089, 0.019102, 0.264006, 0.04071, 0.212324, 0.269086, 0.104357,
                 0.088573, 1, 0.00106, 0.00101, 0.02417, 0.06592, 0.07854, 0.11232,
                 0.71698, 9e-04, 0.9991, 0.00185, 0.00251, 0.8845, 0.11114, 1,
                 0.272, 0.122, 0.238, 0.083, 0.172, 0.057, 0.056, 0.0307, 0.1499,
                 0.1124, 0.1382, 0.0738, 0.2675, 0.2275, 0.4781, 0.5219, 0.002,
                 0.0218, 0.148, 0.2047, 0.1565, 0.2484, 0.2186, 1, 0.00056, 0.00095,
                 0.02329, 0.18889, 0.25475, 0.24896, 0.2826, 1, 0.00139, 0.01601,
                 0.33503, 0.22869, 0.26978, 0.1491, 1, 0.0013, 0.0304, 0.1428,
                 0.2183, 0.1613, 0.3183, 0.1276, 0.9741, 0.0259, 0.0016, 0.0526,
                 0.186, 0.2728, 0.1362, 0.3508, 0.00012, 0.99988, 0.0012, 0.265,
                 0.1431, 0.3064, 0.2843, 0.374, 0.626, 2e-04, 0.0159, 0.0196,
                 0.1324, 0.1615, 0.2626, 0.4078, 0.373, 0.627, 0.00014, 0.00782,
                 0.32967, 0.33832, 0.25242, 0.07163, 1, 0.0015, 0.0997, 0.1687,
                 0.231, 0.1318, 0.2986, 0.0687, 0.2952, 0.7048, 0.014, 0.241,
                 0.221, 0.524, 1, 1, 1, 5.4e-05, 0.007204, 0.992742, 1, 1, 1,
                 1, 1, 1, 1, 1, 0.0075, 0.0425, 0, 1e-04, 0.7576), ratioC = c(6L,
                                                                              6L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 4L, 3L, 3L, 3L, 6L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 3L, 3L, 3L, 3L,
                                                                              3L, 2L, 2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              2L, 2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                                                              0L, 0L, 0L, 0L, 0L, 0L, 4L, 0L, 6L, 2L, 3L, 3L, 0L, 6L, 3L, 3L,
                                                                              3L, 3L, 2L)), row.names = c(NA, -308L), class = "data.frame")



#############################################################################
#############################################################################

## NOT WORKING YET ##
#Add the ISRS formula to the CP_allions table if they exist

addISRS <- function(ISRS_input, CP_allions, threshold) {

    ISRS <- unlist(str_split(ISRS_input, "\n"))

        ISRS_data <- list()

        ISRS <- tibble::tibble(text = ISRS) %>%
            tidyr::separate(text, into = c("Compound_Class", "Adduct_Formula", "Charge"), sep = " ") |>
            dplyr::mutate(Charge = dplyr::case_when(
                Charge == "+" ~ 1,
                Charge == "-" ~ -1))

        for (i in seq_along(ISRS$Adduct_Formula)) {
            Compound_Class <- ISRS$Compound_Class[i]
            Adduct_Formula <- ISRS$Adduct_Formula[i]
            Charge <- ISRS$Charge[i]

            dat <- enviPat::isopattern(isotopes = isotopes,
                                       chemforms = Adduct_Formula,
                                       threshold = threshold,
                                       emass = 0.00054857990924,
                                       plotit = FALSE,
                                       charge = Charge)

            dat <- tibble::as_tibble(dat[[1]], .name_repair = "unique")

            names(dat) <- make.unique(names(dat))

            dat <- dat |>
                dplyr::mutate(abundance = round(abundance, 1)) |>
                dplyr::mutate(`m/z` = round(`m/z`, 6)) |>
                dplyr::mutate(Isotope_Formula = NA) |>
                dplyr::mutate(Molecule_Formula = Adduct_Formula) |>
                dplyr::mutate(Compound_Class = Compound_Class) |>
                dplyr::mutate(Halo_perc = NA) |>
                dplyr::mutate(Adduct_Formula =  Adduct_Formula) |>
                dplyr::mutate(Charge = Charge) |>
                dplyr::mutate(Isotopologue = NA) |>
                dplyr::mutate(Adduct = NA) |>
                dplyr::rename(Rel_ab = abundance) |>
                dplyr::select(Molecule_Formula, Compound_Class, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab)

            ISRS_data[[i]] <- dat

        }
        # combine all elements in list list to get dataframe
        ISRS_data <- do.call(rbind, ISRS_data)
        df <- dplyr::bind_rows(CP_allions, ISRS_data)

    #}
    return(df)
}

