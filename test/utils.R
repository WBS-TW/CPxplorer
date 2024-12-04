# Various utilities and helper functions

create_formula <- function(C, H, Cl, Br, S, O) {
        formula <- paste0(
                case_when(C < 1 ~ paste0(""),
                          C == 1 ~ paste0("C"),
                          C > 1  ~ paste0("C", C)),
                case_when(H < 1 ~ paste0(""),
                          H == 1 ~ paste0("H"),
                          H > 1  ~ paste0("H", H)),
                case_when(Cl < 1 ~ paste0(""),
                          Cl == 1 ~ paste0("Cl"),
                          Cl > 1  ~ paste0("Cl", Cl)),
                case_when(Br < 1 ~ paste0(""),
                          Br == 1 ~ paste0("Br"),
                          Br > 1  ~ paste0("Br", Br)),
                case_when(S < 1 ~ paste0(""),
                          S == 1 ~ paste0("S"),
                          S > 1  ~ paste0("S", S)),
                case_when(O < 1 ~ paste0(""),
                          O == 1 ~ paste0("O"),
                          O > 1  ~ paste0("O", O))
        )
        # Remove any leading or trailing spaces
        str_trim(formula)
}


create_elements <- function(data) {
        # String vector
        string_vector <- c("m/z", "abundance", "12C", "13C", "1H", "2H","35Cl", "37Cl", "79Br", "81Br", "16O", "17O", "18O", "32S", "33S", "34S", "36S")
        
        # Identify columns in the string vector that are not in the data frame
        new_columns <- setdiff(string_vector, names(data))
        
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
        str_trim(formula_iso)
}

# calculate_haloperc <- function(Molecule_Formula) {
#         molform <- rcdk::get.formula(Molecule_Formula)
#         mw <- molform@mass
#         mw_halo <- as_tibble(molform@isotopes) |> 
#                 mutate(mass = as.double(mass)) |> 
#                 mutate(number = as.numeric(number)) |> 
#                 mutate(Halogen = case_when(isoto == "Cl" ~ TRUE,
#                                            isoto == "Br" ~ TRUE,
#                                            isoto == "F" ~ TRUE,
#                                            isoto == "I" ~ TRUE,
#                                            .default = FALSE)) |> 
#                 
#                 filter(Halogen == TRUE) |> 
#                 summarise(mw_halogens = sum(number * mass))  
#         Molecule_Halo_perc <- round(mw_halo$mw_halogens/mw*100, 0)
#         return(Molecule_Halo_perc)
# }

calculate_haloperc <- function(Molecule_Formula) {
        # Regular expression to extract atoms and their counts
        pattern <- "([A-Z][a-z]*)(\\d*)"
        mwtable <- tibble(Atom = c("H", "C", "O", "S", "Cl", "Br"), MW = c(1.00794, 12.011, 15.9994, 32.066, 35.4527, 79.904))
        
        # Extract matches
        matches <- str_match_all(Molecule_Formula, pattern)[[1]]
        
        # Convert to a data frame for clarity
        result <- data.frame(
                Atom = matches[, 2],                     # Element symbols
                Count = as.numeric(matches[, 3])         # Element counts
        )
        
        # Replace missing counts (e.g., implicit "1") with 1
        result$Count[is.na(result$Count)] <- 1
        
        result <- result |> 
                left_join(mwtable, by = "Atom") |> 
                mutate(MW_atoms = MW*Count) |> 
                mutate(Halogen = case_when(Atom == "Cl" ~ TRUE,
                                           Atom == "Br" ~ TRUE,
                                           Atom == "F" ~ TRUE,
                                           Atom == "I" ~ TRUE,
                                           .default = FALSE))
        
        mw <- sum(result$MW_atoms)
        
        mw_halo <- result |> 
                filter(Halogen == TRUE) |> 
                summarise(sum(MW_atoms)) |> 
                as.double()
        
        Molecule_Halo_perc <- round(mw_halo/mw*100, 0)
        return(Molecule_Halo_perc)
}

############################################################################
######### This function generates input for the Envipat function ###########
############################################################################

generateInput_Envipat_normal <- function(data = data, group = group, adduct_ions = adduct_ions, fragment_ions = fragment_ions) {
        
        
        data <- data |> 
                mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                             group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                mutate(Adduct = adduct_ions) |> 
                mutate(Cl = case_when(
                        fragment_ions == "-Cl" ~ Cl-1,
                        fragment_ions == "-HCl" ~ Cl-1,
                        fragment_ions == "+Cl" ~ Cl+1,
                        fragment_ions == "-Cl-HCl" ~ Cl-2,
                        fragment_ions == "-Cl-2HCl" ~ Cl-3,
                        fragment_ions == "-Cl-3HCl" ~ Cl-4,
                        fragment_ions == "-Cl-4HCl" ~ Cl-5,
                        fragment_ions == "-2Cl-HCl" ~ Cl-3,
                        .default = Cl)) |> 
                mutate(H = case_when(
                        fragment_ions == "-H" ~ H-1,
                        fragment_ions == "-HCl" ~ H-1,
                        fragment_ions == "-Cl-HCl" ~ H-1,
                        fragment_ions == "-Cl-2HCl" ~ H-2,
                        fragment_ions == "-Cl-3HCl" ~ H-3,
                        fragment_ions == "-Cl-4HCl" ~ H-4,
                        fragment_ions == "-2Cl-HCl" ~ H-1,
                        .default = H)) |> 
                mutate(Br = case_when(
                        fragment_ions == "+Br" ~ 1,
                        TRUE ~0
                )) |> 
                mutate(Adduct_Formula = case_when(
                        fragment_ions != "+Br" ~ paste0("C", C, "H", H, "Cl", Cl),
                        fragment_ions == "+Br" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br))) |> 
                select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, C, H, Cl)
        
        return(data)
}





# For bromochloroalkanes
generateInput_Envipat_BCA <- function(data = data, group = group, adduct_ions = adduct_ions, fragment_ions = fragment_ions) {
        
        
        data <- data |> 
                mutate(Halo_perc = round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0)) |> 
                mutate(Adduct = adduct_ions) |> 
                mutate(Cl = case_when(
                        fragment_ions == "-Cl" ~ Cl-1,
                        fragment_ions == "-HCl" ~ Cl-1,
                        fragment_ions == "+Cl" ~ Cl+1,
                        fragment_ions == "-Cl-HCl" ~ Cl-2,
                        fragment_ions == "-Cl-2HCl" ~ Cl-3,
                        fragment_ions == "-Cl-3HCl" ~ Cl-4,
                        fragment_ions == "-Cl-4HCl" ~ Cl-5,
                        fragment_ions == "-2Cl-HCl" ~ Cl-3,
                        .default = Cl)) |> 
                mutate(H = case_when(
                        fragment_ions == "-H" ~ H-1,
                        fragment_ions == "-HCl" ~ H-1,
                        fragment_ions == "-Cl-HCl" ~ H-1,
                        fragment_ions == "-Cl-2HCl" ~ H-2,
                        fragment_ions == "-Cl-3HCl" ~ H-3,
                        fragment_ions == "-Cl-4HCl" ~ H-4,
                        fragment_ions == "-2Cl-HCl" ~ H-1,
                        .default = H)) |> 
                mutate(Adduct_Formula = paste0("C", C, "H", H, "Cl", Cl, "Br", Br)) |> 
                select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, C, H, Cl, Br)
        
        return(data)
}



generateInput_Envipat_advanced <- function(data = data, Compounds = Compounds, Adduct_Ion = Adduct_Ion, 
                                  TP = TP, Charge = Charge) {
        
        
        data <- data |> 
                mutate(Adduct_Ion = Adduct_Ion) |> 
                mutate(Charge = Charge) |> 
                mutate(Adduct_Annotation = case_when(
                        TP == "None" ~ paste0("[", Compounds, Adduct_Ion, "]", Charge),
                        .default = paste0("[", Compounds, TP, Adduct_Ion, "]", Charge))) |> 
                mutate(Adduct_Annotation = str_replace(Adduct_Annotation, "\\d$", "")) |> 
                mutate(Compound_Class = Compounds) |> 
                mutate(TP = TP) |> 
                mutate(Cl = case_when(
                        Adduct_Ion == "-Cl" ~ Cl-1,
                        Adduct_Ion == "-HCl" ~ Cl-1,
                        Adduct_Ion == "+Cl" ~ Cl+1,
                        Adduct_Ion == "-Cl-HCl" ~ Cl-2,
                        Adduct_Ion == "-Cl-2HCl" ~ Cl-3,
                        Adduct_Ion == "-Cl-3HCl" ~ Cl-4,
                        Adduct_Ion == "-Cl-4HCl" ~ Cl-5,
                        Adduct_Ion == "-2Cl-HCl" ~ Cl-3,
                        .default = Cl)) |> 
                mutate(H = case_when(
                        Adduct_Ion == "-H" ~ H-1,
                        Adduct_Ion == "-HCl" ~ H-1,
                        Adduct_Ion == "-Cl-HCl" ~ H-1,
                        Adduct_Ion == "-Cl-2HCl" ~ H-2,
                        Adduct_Ion == "-Cl-3HCl" ~ H-3,
                        Adduct_Ion == "-Cl-4HCl" ~ H-4,
                        Adduct_Ion == "-2Cl-HCl" ~ H-1,
                        .default = H)) |> 
                mutate(Br = ifelse(Compound_Class == "BCA", Br, 0)) |> 
                mutate(Br = ifelse(Adduct_Ion == "+Br", Br+1, Br)) |> 
                mutate(O = case_when(
                        TP == "+OH" ~ 1,
                        TP == "+2OH" ~ 2,
                        TP == "+SO4" ~ 4,
                        TP == "-Cl+OH" ~ 1,
                        TP == "-2Cl+2OH" ~ 2,
                        .default = 0
                )) |> 
                mutate(S = case_when(
                        TP == "+SO4" ~ 1,
                        .default = 0)) |> 
                mutate(Adduct_Formula = create_formula(C, H, Cl, Br, S, O))|> 
                rowwise() %>%
                mutate(Molecule_Halo_perc = calculate_haloperc(Molecule_Formula)) |> 
                ungroup() |> 
                select(Molecule_Formula, Molecule_Halo_perc, Charge, Compound_Class, TP, Adduct_Ion, Adduct_Annotation, Adduct_Formula, C, H, Cl, Br, S, O)
        
        return(data)
}
