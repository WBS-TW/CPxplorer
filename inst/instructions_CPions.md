
# Introduction
This app generates a list of mass-over-charges (m/z) for PCAs and structural analogues to investigate potential overlapping m/z during mass spectrometric analysis.
It utilizes the isopattern function of the R package Envipat (Loos et al, Analytical Chemistry, 2015, 87, 5738-5744) to generate isotopic patterns from chemical formula and adducts. The app will also generate a list of m/z for selected adducts which can be used as a transition list in the Skyline software for quantification purposes.
Several adducts have been included and more can be added upon request.  
  
## Instructions    
Choose the parameters in the _Normal settings_ or _Advanced settings_ tab. Press submit and wait for calculation to finish. A table will then be generated with all ions that conform with the initial setting parameters. The table can be exported to excel by clicking on the "Excel" button at the top.  
The _Interfering ions_ tab can be used to check for ions that interfer with each other at the set resolution of the mass spectrometer. Default is set to R=60,000. The plots and tables are interactive and the user can filter the _"Interference at MS res?"_ by clicking on _FALSE_ on the plot legend (and thereby keeping all _TRUE_ ions, which will remove all ions that can be resolved by the set MS resolution).  
The generated ion tables can be exported using the Skyline tab for data processing of MS data.  


## Normal settings tab
  
__C atoms min__ and __C atoms max__: the range of carbon atoms.
  
__Cl atoms min__ and __Cl atoms max__: the range of chlorine atoms.  

__Br atoms min__ and __Br atoms max__: the range of bromine atoms. This is only used if BCA (bromo-chloro alkanes) is choosen as adducts (otherwise this parameter can be ignored).
  
__Add adducts/fragments__: refers to the formula of adducts and/or fragments to generate from a set list of available options. Multiple selections are possible.  
  
__[PCA]__ and __[PCO]__: refers to the main groups of polychlorinated alkanes [PCA] or polychlorinated alkenes (mono olefins) [PCO]. 

__[BCA]__: bromo-chloro alkanes.
  
__[xxx-yy]__ or __[xxx-yy-zz]__: where _xxx_ refers to the main groups, _-yy_ refers to the adduct/fragment ions or _-yy-zz_ which are the consecutive loss fragments. Currently, a limited selection is available but more can be added upon request.  
  
__[ ]-__ and __[ ]+__ refers to the charge of the ion (limited to single charged species, +1 or -1).  
  
_Note that [M+Cl-HCl]- can also be written as [M-H]-_  
  
__Isotope rel ab threshold (%)__: is the threshold for relative abundance for isotopologues for each chemical formula of the adduct/fragment ion. Ions below this threshold will not be included into the generated ion table.
  
__Optional: add ion formula for IS/RS__: input the ion formula for IS or RS if needed. Enter one formula per line. Indicate IS or RS isotopic formula and charge separated by space. _The ion formula should be for the adduct ion_ and not the neutral chemical formula.  
The input for the ion formula is written in three parts, each separated by a `blank space`, the IS or RS, the ion formula where heavy isotopes are indicated by [13]C or [2]H, and the charge (+ or -).  
Example:  
IS [13]C10H16Cl7 -  
RS [13]C12H18Cl9 -  
_(make sure there is no empty last line)_  

The first line will produce the m/z for the [M+Cl]- adduct ion of the IS with the formula ^13^C~10~H~16~Cl~6~, while the second line is for the [M+Cl]- adduct ion for the ^13^C~12~H~18~Cl~8~ RS.  

  
## Advanced settings tab  
Mostly same initial parameters as Normal settings. In advanced settings, there is more flexibility to combine and mix the `Compound Class`, `Adduct`, `Charge`, and `Transformation product`.  
  
### Output table  
  
__Parent_Formula__: the chemical formula of the molecular ion.  

__Halo_perc__: for PCA, PCO chlorinated groups, the chlorination degree (mass percentage of Cl). If BCA, then it is the combined Cl+Br percentage.
  
__Charge__: The charge of the ion.  
  
__Fragment__: The fragment and isotopic type of the ion species.  
  
__Adduct_Formula__: the chemical formula of the adduct/fragment ion.  
  
__Isotopologue__: the isotopologue in relation to the monoisotopic ion.  
  
__Isotope_Formula__: the exact isotopic formula of the adduct/fragment ion.  
  
__m/z__: the mass-over-charge of the adduct/fragment ion.  
  
__Rel_ab__: the relative abundance of the different isotopologues of each adduct/fragment ion.  
  
__12C__, __13C__, __1H__, __2H__, __35Cl__, __37Cl__, __79Br__, __81Br__: the number of atoms for each element
  
## Interfering ions tab  
  
__difflag__, __difflead__: internal calculations for the difference in m/z values between the two nearest ions. If "interference at MS res?" filter has been used, then the previous/next ions might not be shown. 
  
__reslag__, __reslead__: internal calculations for the MS resolution needed to separate the two nearest ions. If "interference at MS res?" filter has been used, then the previous/next ions might not be shown.  
  
__interference__: indicate whether or not the m/z two nearest ions can interfere with each other at the set MS resolution value. _"false"_ means no interference and _"true"_ means there is interference (and therefore the MS resolution cannot resolve these peaks).  


  
  
### Plot outputs  
  
Note: Some bars can exceed 100% in relative abundance in the y-axis, and this indicates that some isotopologues have the exact same ion formula. You can hover over the different segments of these bars to check the overlapping mass ions.

__Hover text__:   
__difflag & difflead (prev and next)__: the difference between the m/z with previous and next ions (axis ordered from lowest to highest m/z).  
__reslag & reslead (prev and next)__: the MS resolution needed to resolve previous and next ion (axis ordered from lowest to highest m/z).  
  
  
## Skyline tab  

This will output a transition list for import in Skyline for integration and quantification. User need to indicate either the _Normal setting_ or _Advanced setting_ as input.  
__Use as Quant Ion__: If _"Most intense"_ then the ion with highest intensity (100%) will be used as the quantification ion. 
The __Interference ions__ calculation MUST be performed before choosing the _"Most intense after interference filtering"_. Ions that have mass interference at the set MS Resolution
will be removed. If no ions are left for a Molecule Name after filtering, then the most intense ions are still kept but [INTERFERENCE] will be added to the Note column.  

__Molecule List Name__: specifies how the compounds are grouped in Skyline. In this case, the naming is by: Compound Class, Carbon Chain Length, Transformation Product (if present). The user can freely change this in the exported spreadsheet for their own needs.  
__Molecule Name__: is the chemical formula of the neutral compound (without adduct/fragment). Beware: Some transformation products can give exactly the same molecular formula and m/z. For example PCA-Cl+OH will be exactly same as PCA-H+OH with one less chlorine atom. 
C10H18Cl4(-Cl+OH) -> C10H19Cl3O and C10H19Cl3(-H+OH) -> C10H19Cl3O.  
__Precursor Charge__: is the chosen ionization mode.  
__Label Type__: Quan is the highest intensity isotopologue of each Compound. Can be changed in the exported spreadsheet.  
__Precursor m/z)__: is the m/z after ionization and adduct/fragment formation.  
__Explicit Retention Time/Window__: used by Skyline to integrate peaks. Can be left empty.  
__Note__: Internal information on ion formula and isotopic abundance. DO NOT change this (used by CPquant).  





  



