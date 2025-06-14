library(dada2)
library(ggplot2)
set.seed(100)
#for Windows
setwd("~/personal analysis/R")
seq_folder <- "~/Rawdata" #change to the directory containing the fastq raw files
analysis_path <- "~/personal analysis/R"

list.files(seq_folder)
file_folders <-list.files(seq_folder)
file_folders <- file_folders[1:5] #remove two files from the list based on the list.files results
my_16S_folder_path <- vector(length(file_folders), mode="character")
fnFs <- vector(length(file_folders), mode="character")
fnRs <- vector(length(file_folders), mode="character")
for (i in 1:length(file_folders)) {
  my_16S_folder_path[i] <- paste(seq_folder, file_folders[i], sep="/")
  fnFs[i] <- list.files(my_16S_folder_path[i], pattern="[0-9A-Za-z]_1.fq.gz", recursive = TRUE, full.names = TRUE)
  fnRs[i] <- list.files(my_16S_folder_path[i], pattern="[0-9A-Za-z]_2.fq.gz", recursive = TRUE, full.names = TRUE)
}


sample.names <- file_folders
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

#use MicrobiomeAnalystR data2_utilities.R for dada2
#run data2_utilities_revised.R first


setParametersRes <- setParameters(file_compressed = TRUE, 
              OS_is_windows = TRUE)

.plotQualityProfileLoop(f_F = fnFs, #forward reads file
                        f_R = fnRs,#reverse reads file
                        sn = sample.names, # sample name
                        plot_format = "pdf", # pdf, tiff, ....
                        fd = analysis_path)

#DON'T use the function seAanityCheck, because raw files are not in the same folder
#use the following script instead
seqSanityCheckRes <- list(fnFs = fnFs,
                          fnRs = fnRs,
                          sn = sample.names)


processRawSeqRes <- processRawSeq(setParametersRes = setParametersRes, # results from setParameters
              seqSanityCheckRes = seqSanityCheckRes, # results from seqSanityCheck
              reads_trim_length_F_R = '', #for V3-V4 no need for this setting #retained reads length for forward and reverse,
              #depending on the quailty of reads, the quailty graph can be obtained by seqSanityCheck results;
              plot_format = "pdf")


constructSeqTabRes <- constructSeqTab(setParametersRes = setParametersRes, # results from setParameters
                                      processRawSeqRes = processRawSeqRes)


#check the dataset URL first in data2_utilities.R
assignTaxRes <- assignTax(constructSeqTabRes = constructSeqTabRes, #results from constructSeqTab
          ref_db = "silva")

constructPhyloTreeRes <- constructPhyloTree(constrcutSeqTabRes = constrcutSeqTabRes)
