
# Introduction    
CPquant uses the deconvolution process proposed by Bogdal et al (Anal Chem, doi/10.1021/ac504444d) to estimate 
the relative composition needed from different standards to match the measured homologue pattern of samples. 
The underlying calculations are based on the CPCrawler R script by Perkons et al (Food Chem, doi/10.1016/j.foodchem.2019.125100).
In CPquant, the deconvolution is performed using the nnls package (https://cran.r-project.org/web/packages/nnls).
  
  
## Input file  
The input excel file should be from the Skyline results table. It must include the following column with the names:  
`Replicate Name`: sample name  
`Sample Type`: Unknown (which is the sample to be quantified), Blank, Standard, Quality Control  
`Molecule List`  
`Molecule`  
`Area`  
`Mass Error PPM`  
`Isotope Label Type`  
`Chromatogram Precursor M/Z`  
`Analyte Concentration`: For standards only. This is the standard concentrations/amounts. This column could be in concentration or weight/amount unit depending on the user input. It will affect the final quantification unit.  
`Batch Name`: For standards only. This will determine which standards that belongs to a concentration series as well as which carbon chain groups to quantify with the standard.  
The naming of the Batch Name should be: CarbonGroups_StandardName. An underscore is separator for the carbon chain group and standard name.  
Examples:  
C10-C13_StandardA. This standard will then be used to quantify carbon chains C10, C12, C13. This belongs to the StandardA which can be in different concentrations for the calibration series.  
C14_52%Cl. This standard will only quantify C14 carbon chains. It specifies 52% chlorine content (this information is not needed for quantification).  
  
  
## Quantification Inputs tab    
After loading the excel, allow for the Area plot to show up before pressing the "Proceed" button, otherwise error will occur.  
After loading the data, the user can choose the options:  
__Subtraction by blank?__: If "Yes, by avg area of blanks", then the area for each Molecule will be subtracted with the average of all blank samples.  
__Correct with RS area?__: If "Yes", then the area of each Molecule will the normalized to the recovery standard (RS) area for each sample.  
__Calculate recovery?__: If "Yes", requires samples with the `Sample Type` designated as "Quality Control" that include the same concentrations of IS and RS.  
__Calculate MDL?__: If "Yes", then calculates the method detection limits based on blank samples.  
If no blank subtraction then MDL = avg + 3 * standard deviation of blank samples.  
If blank subtraction then MDL = 3 * standard deviation of blank samples.  
__Types of standards__: Currently only have option to use mixtures and single chain standards to perform deconvolution. More option can be added later for other quantification strategies.  

  
  
__Remove samples from quantification?__: You can select samples to be removed before quantification process.  
__Keep the the calibration curves above this rsquared__: Keeping all calibration curves for every homologue groups for each standard that are above this R2 value. 
This will remove all homologue groups that do not show linearity within the standard calibration concentrations, thus remove their contribution to the deconvolution. 
Default is 0.8 but can be changed accordingly by the user.  
  
CLICK ON __Proceed__ will quantify the samples based on the deconvolution process and the results will show in the different tabs.  
  
_Quantification process_: The process starts by creating calibration curves for each carbon chain group for each standard mixture.
The Batch Name in the excel file determines which carbon chain group to be included for each standard mixture. A linear regression will be fitted and the slope is used as the response factor (RF).
If the R-squared of the fit for a homologue group for a standard series (calibration curve) is below the user input threshold (modified in the first tab), then the homologue group in that standard is not considered for subsequent quantification.  



## Input summary    
### Choose tab  
Might take some time before results show up here so be patient.  
__Std Calibration Curves__: The calibration curves for different standards will be shown. Only those with rsquared above the initial cutoff will be shown.  
__Quan to Qual ratio__: Violin plots showing the ratio Quan/Qual area to detect outliers and thus help in assessing quality of data.  
  
  
## Quantification summary  
__Export all results to Excel__: export all results from the quantification to an excel file with different sheets.  
__Quantification table__: this shows the quantification results directly in a table. The unit of the quantification depends on the design of concentration or weight amount specified by the user.  
__Contributions from standards to deconvoluted homologue pattern__: this plot shows how much each standard contributes to the reconstructed homologue group pattern.  
  
  
## Homologue Group Patterns  
  
Plots the relative distribution (relative area) of the samples.  
__All Samples Overview__: gives a quick overview on homologue group patterns of all samples in a static plot.  
__Samples Overlay__: overlays all selected samples in one plot.
__Samples Panels__: plots one panel for each selected sample. Also compares the relative distribution of homologue groups of the sample with the reconstructed pattern 
by the deconvolution process (as scatter lines in the Deconvoluted Distribution legend group).
BE AWARE: CURRENTLY THE COLORS OF THE CARBON CHAIN GROUPS DOES NOT MATCH BETWEEN DIFFERENT SAMPLES  
  
  
## QA/QC  
Various QA/QC results will show depending on the choices in the input tab. These include recovery and MDL calculations.  








