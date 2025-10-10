
# CPquant
  
  
## Introduction    
CPquant uses the deconvolution process proposed by Bogdal et al (Anal Chem, doi/10.1021/ac504444d) to estimate 
the relative composition needed from different standards to match the measured homologue pattern of samples. 
The underlying calculations are based on the CPCrawler R script by Perkons et al (Food Chem, doi/10.1016/j.foodchem.2019.125100).
In CPquant, the deconvolution is performed using the nnls package (https://cran.r-project.org/web/packages/nnls).  

- CPquant quantification works with both single chain and mixture standards.  
- We recommend to use 5 calibration levels for each standard. These needs to be named in the `Batch Name` column and their concentration levels `Analyte Concentration` added according to below instructions from the input file which is exported from Skyline.  
- If recovery needs to be calculated, then a _Quality Control_ sample needs to be added.

__The calculated concentrations are for those in the extract. The user can then export the results to excel and perform additional calculations to derive the concentrations in the samples__
  
  
## Input file  
The input excel file should be exported from the Skyline results table. It should include the following column with the names:  
`Replicate Name`: sample name  
`Sample Type`: the following characters can be used, _Unknown_ (which is the sample to be quantified), _Blank_ (field blanks and procedural blanks are not distinguished), _Standard_ (standard used for quantification), _Quality Control_ (standard/sample used to determine the recovery).  
`Batch Name`: For standards only (leave blank for Unknown). This will determine which standards that belongs to a calibration series as well as which carbon chain groups 
to quantify with the standard.  
The naming of the Batch Name should be: CarbonGroups_StandardName. An underscore is a separator for the carbon chain group and standard name.  
Example A: C10-C13_StandardA. This standard will then be used to quantify carbon chains C10, C12, C13 (the hyphen specify the range of carbon chains). 
This belongs to the StandardA which can be at different Analyte Concentration for the same calibration series.  
Example B: C14_52%Cl. This standard will only quantify C14 carbon chains. It specifies 52% chlorine content (although this information is not needed for quantification).  
`Molecule List`: compounds used internal standards are denoted `IS`, recovery standards as `RS` (also called volumetric standard).  
`Molecule`: PCA homologue group.  
`Area`: integrated area from Skyline.  
`Analyte Concentration`: For standards only. This is the standard concentrations/amounts. 
This column could be in concentration or weight/amount unit depending on the user input. It will affect the final quantification unit.  
`Mass Error PPM`: might be exported from Skyline but currently not used by CPquant.  
`Isotope Label Type`: Quan or Qual.  
`Chromatogram Precursor M/Z`: the m/z values of the ion (not used by CPquant). 
`Sample Dilution Factor`: indicate the dilution (>1) or concentration factor (<1). Default from Skyline is 1 (no dilution).  
`Transition Note`: internal information transferred from CPions used for calculating the correct isotope ratio.  
  
## Quantification Inputs tab  
__Import excel file from Skyline__: This is the excel file from the Report export function of Skyline.  
__Concentration unit__: an optional input to indicate the concentration (or amount) unit of the Analyte Concentration (e.g. ng/mL or ng).   
After loading the excel, allow for the Area plot to show up before pressing the "Proceed" button, otherwise error will occur.  

After loading the data, the user can choose the options:  
__Choose ion for quantification__: "Quan only" only uses the signal from quantification ion, and "Sum Quan+Qual" use the sum of the Quan and all Qual ions for quantification.  
__Subtraction by blank?__: If "Yes, by avg area of blanks", then the area for each Molecule will be subtracted with the average area of all blank samples.  
__Correct with RS area?__: If "Yes", then the area of each Molecule will the normalized to the recovery standard (RS) area for each sample.  
__Calculate recovery?__: If "Yes", requires samples with the `Sample Type` designated as "Quality Control" that include the spiked concentrations of 
IS and RS corresponding amounts/concentrations.  
__Calculate MDL?__: If "Yes", then calculates the method detection limits based on blank samples.  
If no blank subtraction then MDL = avg + 3 * standard deviation of blank samples.  
If blank subtraction then MDL = 3 * standard deviation of blank samples.  
__Types of standards__: Currently only have option to use mixtures and single chain standards to perform deconvolution. More option can be added later for other quantification strategies.  
  
  
__Remove samples from quantification?__: select samples to be removed before quantification process.  
__Keep the the calibration curves above this rsquared__: remove calibration curves for every homologue group in each standard below this R2 value (goodness of calibration fit). 
This will remove all homologue groups that do not show linearity within the standard calibration levels, thus remove their contribution to the deconvolution. 
Default is 0.8 but can be changed accordingly by the user.  
  
__Proceed__: pressing this button will quantify the samples based on the deconvolution process and the results will show up in the different tabs.  
  

### Quantification process  
The process starts by creating calibration curves for each carbon chain group for each standard mixture. 
The Batch Name in the excel file determines which carbon chain group to be included for each standard mixture. A linear regression will be fitted and the slope is used as the response factor (RF).
If the R-squared of the goodness of calibration fit for a homologue group for a standard series (calibration curve) is below the user input threshold (modified in the first tab), then the homologue group in that standard is not considered for subsequent quantification.  



## Input summary    
### Choose tab  
The display might take some time before results show up here so be patient.  
__Std Calibration Curves__: The calibration curves for different standards will be shown. Only those with rsquared above the initial cutoff will be shown.  
__Removed from Calibration__: A table showing individual homologue groups from specific standards that are removed from the quantification process, due to negative RF or calibration curve R2 values below limit.   
__Quan to Qual ratio__: Violin plots showing the ratio Quan/Qual area to detect outliers and thus help in assessing quality of data.  
__Measured vs Theor Quan/Qual ratio__: Plot showing the measured Quan/Qual ratio divided by the theoretical Quan/Qual ratio. Ideally, the ratio should be around 1. Outlier ratios (<0.3 or >3) are marked in red.  
  
  
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
  








