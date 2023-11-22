Info:
*****
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: READ_ME

get started
***********

### requirements raw files ###
# do not change names of the raw files, structured name eg. 180616s01 
# filename must start with number, no spaces, no -,#,... symbols
# Importing data into R pipeline
- use samples after conditoning step
- paste all biological samples (samples, QCs, optional: STDs/ISTDs) in folder named 'bio'
- paste all blanks in folder named 'blank' (solvent samples eg. MeOH, ACN)
  ! min. 1 blank sample is needed/used for pre-processing with LC-MS, not needed for direct MS
- optional: zip folder bio
- optional: zip folder blank
- paste (zipped) folder(s) in 'Pipeline_metabolomics/Data/Input' folder for analysis
  ! LC-MS: must be folders blank and bio in Input directory for pre-processing, containing at least 1 sample in each folder 
  ! MS: min presence of folder bio, containing at least 1 sample


################################
### file conversion MS Waters via Rstudio (windows, outside R pipeline)###
# => REIMS file converter microscript in Rstudio
# this microscript runs outside the R pipeline 
# requirments: Windows OS and waters folder present on computer
# to run, open in Rstudio, adjust the 'adjustments' below if needed and click "Source" to run the whole script
# no progressbars are shown, as long as red running symbol is visible, the script is calculating
# check file conversion succesful after finished (no symbols present)


### file conversion LC-MS via GUI (optimalisation)###
# A) convert .raw files from (Q-)Exacutive into mzXML files
#! check NO spaces in filenames!
- open MSConvertGUI on windows
- select folder saved .raw files
- press 'add' button
- choose mzXML as outfut format
- remove filter present
- add filter: peakpicking > vendor with ms levels from 1
- add filter: subset > scan polarity > pos/neg
- leave other parameters on default (64-bit, use zip, ...)
- press 'start' button



### Adjust settings for your analysis ###
"configuration file" = all settings for the code to run.
- copy Configuration(reserve, see input).R into you input
- rename to Configuration.R
- Open Configuration.R 
- adjust settings under 'Configuration > Adjustments' according to the options given under 'Configuration > Options'
- projectname is structured and as short as possible: YYYYMMDD_initials_description_ionisationMode
  eg. 20191031_MDG_exp1_pos
  seperate with "_", don't use symbols (like +,-), spaces, ... 
  description part is optional, max 1 word
- polarity is needed for the pre-processing modules, options are "positive", "negative"
- user_comment is an optional, short comment that will be displayed on the first page of report; eg  explain (Multiple)Comparisons
- when only part of the pipeline is run: adjust the global settings (eg. projectname) and the selected parts
- Save Configuration.R
- close Configuration.R
