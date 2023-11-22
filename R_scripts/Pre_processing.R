# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Part I: pre-processing


##########R Pipeline - Part I: pre-processing##########
print(Sys.time())
start_time <- Sys.time()
print("R pipeline - Part I: pre-processing - start!")
# Part I: pre-processing

## data_loading
#read all mzXML files in folder:
#https://stackoverflow.com/questions/9564489/opening-all-files-in-a-folder-and-applying-a-function
setwd(path_data_in)
filenames <- list.files(path_data_in, pattern="*.mzXML", recursive = TRUE)

# unzip files if needed
path_data_in_bio <- file.path(path_data_in, 'bio.zip') 
path_data_in_blank <- file.path(path_data_in, 'blank.zip')
if(length(filenames) == 0){
  unzip(path_data_in_bio, exdir=path_data_in)
  unzip(path_data_in_blank, exdir=path_data_in)
  filenames <- list.files(path_data_in, pattern="*.mzXML", recursive = TRUE)
}



## Set peak picking parameters
if(PARAMETERS_PREPROCESSING == OPTIMIZED_PEAK_PICKING){ 
  ## Optimisation of peak picking parameters by using natural, stable 13C isotopes
  # optimized parameters for peak picking, retention time correction and peaks filling.
  #https://bookdown.org/yufree/Metabolomics/demo.html
  #https://www.bioconductor.org/packages/devel/bioc/vignettes/IPO/inst/doc/IPO.html
  library(IPO)
  peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave') #for high resolution MS
  peakpickingParameters$ppm <- 6 #same as MZstep parameter
  peakpickingParameters$min_peakwidth <- c(5,20) #45s gives error
  peakpickingParameters$noise <- 15000 #moet >0 want anders superlang berekenen... neem zelfde als pdf
  
  resultPeakpicking <- optimizeXcmsSet(files = filenames[c(1,2,3)], params = peakpickingParameters, nSlaves = 1, subdir = NULL, plot=F) #optimize using first 3 samples
  optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset
  
  retcorGroupParameters <- getDefaultRetGroupStartingParams()
  resultRetcorGroup <- optimizeRetGroup(xset = optimizedXcmsSetObject, params = retcorGroupParameters, subdir = NULL)
  
  #writeRScript(resultPeakpicking$best_settings$parameters, resultRetcorGroup$best_settings)
  para <- capture.output(writeRScript(resultPeakpicking$best_settings$parameters, resultRetcorGroup$best_settings), type = "message")
  
  setwd(path_data_out)
  save(para,file = paste(name_project, '_para.RData', sep=""))
  write_dataframe_as_txt_file(para, paste(name_project, '_Settings_pre-processing.txt', sep=""))
  dev.off() #only needed if perform optimize param
  setwd(path_data_in)
  
  # settings
  library(stringr)
  peakwidth <- as.numeric(unlist(str_extract_all(para[grepl('peakwidth',para)],'\\d+\\.*\\d*'))) #expected approximate peak width in chromatographic space. Given as a range (min, max) in seconds.
  ppm <- as.numeric(unlist(str_extract_all(para[grepl('ppm',para)],'\\d+'))) #maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition.
  noise <- as.numeric(unlist(str_extract_all(para[grepl('noise',para)],'\\d+'))) #minimum intensity required for centroids to be considered in the first analysis step 
  snthresh <- as.numeric(unlist(str_extract_all(para[grepl('snthresh',para)],'\\d+'))) #signal to noise ratio cutoff.
  mzdiff <- as.numeric(unlist(str_extract_all(para[grepl('mzdiff',para)],'\\d+\\.*\\d*'))) #minimum difference in m/z dimension for peaks with overlapping retention times; can be negatove to allow overlap
  prefilter <- as.numeric(unlist(str_extract_all(para[grepl('prefilter',para)],'\\d+\\.*\\d*'))) #c(k, I), Mass traces are only retained if they contain at least k peaks with intensity >= I
  integrate <- as.numeric(unlist(str_extract_all(para[grepl('integrate',para)],'\\d+'))) #integrate = 1 peak limits are found through descent on the mexican hat filtered data
  profStep <- round(as.numeric(unlist(str_extract_all(para[grepl('profStep',para)],'\\d+\\.*\\d*'))),1)
  center <- as.numeric(unlist(str_extract_all(para[grepl('center',para)],'\\d+')))
  response <- as.numeric(unlist(str_extract_all(para[grepl('response',para)],'\\d+')))
  gapInit <- as.numeric(unlist(str_extract_all(para[grepl('gapInit',para)],'\\d+\\.*\\d*')))
  gapExtend <- as.numeric(unlist(str_extract_all(para[grepl('gapExtend',para)],'\\d+\\.*\\d*')))
  factorDiag <- as.numeric(unlist(str_extract_all(para[grepl('factorDiag',para)],'\\d+')))
  factorGap <- as.numeric(unlist(str_extract_all(para[grepl('factorGap',para)],'\\d+')))
  localAlignment <- as.numeric(unlist(str_extract_all(para[grepl('localAlignment',para)],'\\d+')))
  bw <- as.numeric(unlist(str_extract_all(para[grepl('bw',para)],'\\d+\\.*\\d*'))) #bandwidth (standard deviation ot the smoothing kernel) 
  mzwid <- as.numeric(unlist(str_extract_all(para[grepl('mzwid',para)],'\\d+\\.*\\d*'))) #width of overlapping m/z slices to use for creating peak density chromatogramsand grouping peaks across samples
  #error because x.xx is correct but not 1e-4, becomes list [1, 04], solution:
  if(length(mzwid) != 1){
    mzwid <- mzdiff
  } 
  minfrac <- as.numeric(unlist(str_extract_all(para[grepl('minfrac',para)],'\\d+\\.*\\d*'))) #minimum fraction of samples in at least one sample group in which the peaks have to be present to be considered as a peak group 
  minsamp <- as.numeric(unlist(str_extract_all(para[grepl('minsamp',para)],'\\d+'))) #minimum number of samples in at least one sample group in which the peaks have to be detected to be considered a peak group (feature).
  max <-  as.numeric(unlist(str_extract_all(para[grepl('max',para)],'\\d+')))
}
if(PARAMETERS_PREPROCESSING == FIXED_PEAK_PICKING){ 
  #default from G4M slides for peak picking, rt align and filling
  ppm <- 6 #nstrum setup QE orbitrap ms: 10ppm
  peakwidth <- c(5,45) #instrum setup QE orbitrap ms: 15sec
  snthresh <- 10
  noise <- 15000 #toek todo eval 10+4, 10+5, 10+6 for best balance speed/sensitivity
  mzdiff <- 0.05
  prefilter <- c(3,1000)
  mzwid <- 0.05
  bw <- 30  # peak grouping
  #instrum setup QE orbitrap ms resolution= 140000
  
  setwd(path_data_out)
  save(ppm, peakwidth, snthresh, noise, mzdiff, prefilter, mzwid, bw ,file = paste(name_project, '_para.RData', sep=""))
  settings <- c(ppm, peakwidth, snthresh, noise, mzdiff, prefilter, mzwid, bw)
  write_dataframe_as_txt_file(settings, paste(name_project, '_Settings_pre-processing.txt', sep=""))
  setwd(path_data_in)
}



## feature detection and filtering 
#https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html
#https://rdrr.io/search?package=xcms&repo=bioc&q=plot #for plots
library(xcms)
setwd(path_data_in) #anders error, needs to be in loc of files
xset <- xcmsSet(files = filenames, ppm=ppm, mzdiff=mzwid, noise=noise, snthresh=snthresh, 
			method="centWave", fitgauss=TRUE, nSlaves=8) #todo optimize param

xraw <- getXcmsRaw(xset)
setwd(path_data_out)
pdf(paste(name_project,'_TIC_before_align.pdf', sep=""),height=5,width=7)
plotTIC(xraw) #Plot chromatogram of total ion count. 
dev.off()
setwd(path_data_in)

## RT alignment (adjuct according to retention time)
setwd(path_data_out)
pdf(paste(name_project,'_rt_correction_obiwarp_deviation.pdf', sep=""),height=5,width=7)
xset2 <- retcor(xset, method="obiwarp", plottype = "deviation")
dev.off()
pdf(paste(name_project,'_rt_correction_obiwarp_deviation(2).pdf', sep=""),height=5,width=7)
plotrt(xset2) #Use corrected retention times for each sample to calculate retention time deviation profiles and plot each on the same graph. 
dev.off()
setwd(path_data_in)

if(CODE_RUN_MODE == CODE_AUTORUN){
  rm(xset)
  gc()  #It can be useful to call gc() after a large object has been removed, as this may prompt R to return memory to the operating system. gc() also return a summary of the occupy memory.
}

xraw2 <- getXcmsRaw(xset2)

## feature grouping 
xset3 <- group(xset2, method="nearest") #todo optimize: , bw = bw, minfrac = 0.2, minsamp = 1, mzwid = mzwid, max = 50, sleep = 0) #todo: also see if density/nearest is better method?
if(CODE_RUN_MODE == CODE_AUTORUN){
  rm(xset2)
  gc()
}

## fill chrom peaks
xset4 <- fillPeaks(xset3)
if(CODE_RUN_MODE == CODE_AUTORUN & PLOTTING != SAVE_PLOTS){
  rm(xset3)
  gc()
}
reporttab <- diffreport(xset4, filebase = paste(path_data_out, "/diffreport_", name_project, sep=""), mzdec=4, eicmax=5000, metlin = FALSE) #todo optimize param
#statisitcs on bio vs blank

## feature annotation
library(CAMERA)
gc()
#Error: `geom` must be either a string or a Geom object, not an S4 object with class xcmsSet (only with carol/dorain, not with beata)
an <- annotate(xset4, pval=0.05, #, nSlaves=8, calcIso=TRUE, calcCaS=FALSE, maxcharge=3, maxiso=4, minfrac=0.5, #todo under this selection only ]- or ]+
               ppm=10, mzabs=0.015, quick=FALSE, psg_list = NULL, rules = NULL, polarity = POLARITY) #todo optimize param

if(CODE_RUN_MODE == CODE_AUTORUN){
  rm(xset4)
  gc()
}
diffreport1 <- getPeaklist(an)
if(CODE_RUN_MODE == CODE_AUTORUN){
  rm(an)
}
diffreport1 <- cbind(diffreport1[,c("isotopes", "adduct", "pcgroup")], reporttab) #intesities at the end of df


## make variableMetadata, add CompID nrs, export 
variableMetadata <- NULL
variableMetadata$CompID <- rownames(diffreport1)
variableMetadata$MZ <- diffreport1$mzmed
variableMetadata$Time <- diffreport1$rtmed
variableMetadata$Time <- unlist(lapply(variableMetadata$Time, function(x) x/60)) #put RT in mins
variableMetadata <- cbind(variableMetadata, diffreport1)


## additional blank and noise removal after PP steps
if(BLANK_NOISE_REMOVAL == BLANK_NOISE_REMOVAL_ON){
  #remove noise (from BLANKA public https://link.springer.com/article/10.1007/s13361-019-02185-8), 
  #calc 5% lowest intensity peaks, calc average noise I, noiseI * S/N thresh, rm all peaks below this intesnit
  #https://github.com/gtluu/blanka/blob/master/blanka_lcms.py
  #do accros all bio+blank samples by using average thereof, also use median instead #not!
  variableMetadata$AverageI <- apply(variableMetadata[,20:ncol(variableMetadata)], 1, function(x) median(x))
  baseline_noise <- quantile(variableMetadata$AverageI, .05) 
  variableMetadata <- variableMetadata[variableMetadata$AverageI >= (baseline_noise * 4), ] #snthresh default 10 is too much! pick 4 like publication, especially since already removed blank>bio and after PP
  variableMetadata <- variableMetadata[1:(length(variableMetadata)-1)]
  
  #remove variables with intensity higher in blanks > bio
  variableMetadata <- variableMetadata[variableMetadata$tstat < 0, ] #keep alle the CompIDs with negative t-statistic ('POS if graeter in 2nd group = blanks)
  #no min fold threshold between bio and blank, only check different using tstat #not good to use
  
  #remove variables present in blank sample (from BLANKA public https://link.springer.com/article/10.1007/s13361-019-02185-8)
  #use last sample == blank for the extract, will rm a lot so use threshold ipv only keep feats with blank 0 value. too hard to not rm valuable low i peaks
  #variableMetadata <- variableMetadata[variableMetadata[,ncol(variableMetadata)] < baseline_noise / 4, ] #not!
  
  #with a QC/blank ratio <3
  #variableMetadata <- variableMetadata[variableMetadata$fold < 0, ] #no, will be none left. too strict. keep only tstat
  
  #rename CompIds from 1->n
  #variableMetadata$CompID <- c(1, nrow(variableMetadata))   #skip: for eic plot evaluation
}


## save variablemetadata
setwd(path_data_in)
write_dataframe_as_txt_file(variableMetadata, paste(name_project, '_variableMetadata.txt', sep="")) #input for R_pipeline part II: statistical analysis
setwd(path_data_out)
write_dataframe_as_txt_file(variableMetadata, paste(name_project, '_variableMetadata_output_pre-processing.txt', sep="")) #copy for user




## make additional plots
if (PLOTTING == SAVE_PLOTS){
  setwd(path_data_out)
  
  #feat detection
  pdf(paste(name_project,'_AIC_before_align.pdf', sep=""),height=5,width=7)
  plotChrom(xraw) #plot averaged or base peak extracted ion chromatograms over a specified mass range.
  dev.off()
  
  pdf(paste(name_project,'_logintensity_mz_rt_before_align.pdf', sep=""),height=5,width=7)
  image(xraw) #Create log intensity false-color image of plotted with m/z and retention time axes 
  dev.off()
  
  pdf(paste(name_project,'_mz_rt_before_align.pdf', sep=""),height=5,width=7)
  plotRaw(xraw) #Produce a scatterplot showing raw data point location in retention time and m/z. 
  dev.off()
  
  pdf(paste(name_project,'_intensity_mz_before_align.pdf', sep=""),height=5,width=7)
  plotSpec(xraw) #plot mass spectra over a specified retention time range. 
  dev.off()
  
  pdf(paste(name_project,'_EIC_before_align.pdf', sep=""),height=5,width=7)
  plotEIC(xraw) #Plot extracted ion chromatogram for m/z values of interest. The raw data is used in contrast to plotChrom which uses data from the profile matrix. 
  dev.off()
  
  #rt align
  pdf(paste(name_project,'_TIC_RTalign.pdf', sep=""),height=5,width=7)
  plotTIC(xraw2) #Plot chromatogram of total ion count. 
  dev.off()
  
  pdf(paste(name_project,'_AIC_RTalign.pdf', sep=""),height=5,width=7)
  plotChrom(xraw2) #plot averaged or base peak extracted ion chromatograms over a specified mass range.
  dev.off()
  
  pdf(paste(name_project,'_logintensity_mz_rt_RTalign.pdf', sep=""),height=5,width=7)
  image(xraw2) #Create log intensity false-color image of plotted with m/z and retention time axes 
  dev.off()
  
  pdf(paste(name_project,'_mz_rt_RTalign.pdf', sep=""),height=5,width=7)
  plotRaw(xraw2) #Produce a scatterplot showing raw data point location in retention time and m/z. 
  dev.off()
  
  pdf(paste(name_project,'_intensity_mz_RTalign.pdf', sep=""),height=5,width=7)
  plotSpec(xraw2) #plot mass spectra over a specified retention time range. 
  dev.off()
  
  pdf(paste(name_project,'_EIC_RTalign.pdf', sep=""),height=5,width=7)
  plotEIC(xraw2) #Plot extracted ion chromatogram for m/z values of interest. The raw data is used in contrast to plotChrom which uses data from the profile matrix. 
  dev.off()
  
  #Use "democracy" to determine the average m/z and RT deviations for a grouped xcmsSet, and dependency on sample or absolute m/z 
  pdf(paste(name_project,'_plotQC1.pdf', sep=""),height=5,width=7)
  plotQC(xset3, what="mzdevhist") #histogram of mz deviations.Should be gaussian shaped. If it is multimodal, then some peaks seem to have a systematically higher m/z deviation
  dev.off()
  
  pdf(paste(name_project,'_plotQC2.pdf', sep=""),height=5,width=7)
  plotQC(xset3, what="rtdevhist") #histogram of RT deviations. Should be gaussian shaped. If it is multimodal, then some peaks seem to have a systematically higher RT deviation
  dev.off()
  
  try({
    pdf(paste(name_project,'_plotQC3.pdf', sep=""),height=5,width=7)
    plotQC(xset3, what="mzdevmass") #Shows whether m/z deviations are absolute m/z dependent, could indicate miscalibration 
    dev.off()
    
    pdf(paste(name_project,'_plotQC4.pdf', sep=""),height=5,width=7)
    plotQC(xset3, what="mzdevtime") #Shows whether m/z deviations are RT dependent, could indicate instrument drift
    dev.off()
  })
  
  pdf(paste(name_project,'_plotQC5.pdf', sep=""),height=5,width=7)
  plotQC(xset3, what="mzdevsample") #median mz deviation for each sample, indicates outliers 
  dev.off()
  
  pdf(paste(name_project,'_plotQC6.pdf', sep=""),height=5,width=7)
  plotQC(xset3, what="rtdevsample") #median RT deviation for each sample, indicates outliers
  dev.off()
  
  setwd(path_data_in)
}
if (PLOTTING == DONT_SAVE_PLOTS){
  print("No plots are printed")
}


print("R pipeline - Part I: pre-processing - done!")
print(Sys.time())
end_time <- Sys.time()
print(end_time - start_time)
#
#####################
