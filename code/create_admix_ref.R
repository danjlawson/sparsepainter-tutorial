#!/usr/bin/env Rscript
## Creates a reference panel for admixture estimation from a leave-one-out chromosome painting (made with SparsePainter or ChromoPainter)

if(!exists("args") || class(args)!="character"){
    args <- commandArgs(trailingOnly = TRUE)
}

if(length(args)<1){
    cat("Usage: Rscript create_admix_ref.R <options> <inputfile> <popidfile> <outputfile>\n")
    cat("Creates a reference panel for admixture estimation from a leave-one-out chromosome painting (made with SparsePainter or ChromoPainter)\n")
    cat("inputfile: A file containing the chromosome painting results, with individuals as rows and donor populations as columns. The first column should be the individual IDs.\n")
    cat("popidfile: A file containing the population IDs for each individual, with two columns: the first column should be the individual IDs (matching those in the inputfile), and the second column should be the population IDs.\n")
    cat("outputfile: The file to write the reference panel to, in a format suitable for admixture estimation (e.g., with nnls_panel.R)\n")
    cat("OPTIONS:\n")
    cat("-v : Verbose mode (default: FALSE)\n")
    cat("-c <n>: Number of times to expect each individual in the input file (default: 1). If >1, looks for form _<i-1> at the end of individual IDs to identify replicates. SparsePainter produces -n 2 format data for diploids.\n\n") 
    cat("Example usage: Rscript create_admix_ref.R -v sparsepaint/test.combined.chunklength.txt testcombined.pop.ids test_admix_panel.txt\n")
    stop("Must provide 1 or more arguments")
}

i=1
tempargs=args
copies=1
verbose=FALSE

while(i<=length(tempargs)) {
  if(tempargs[i]=="-v"){
        verbose=TRUE
        tempargs=tempargs[-i]
    } else if (tempargs[i]=="-c"){
        i=i+1
        copies=as.numeric(tempargs[i])
        tempargs=tempargs[-c(i-1,i)]
    } else {
        i=i+1
    }
}
files = tempargs
if(length(files)!=3){
    stop("Must provide exactly 3 file names: input, popidfile, and output")
}

source("../code/FinestructureLibrary.R") # read in the R functions, which also calls the needed packages

data <- as.matrix(read.table(files[1], header=TRUE, row.names=1))
if(verbose) print(paste("Read data with",nrow(data),"individuals and",ncol(data),"donor populations"))
popids <- read.table(files[2], header=FALSE, row.names=1)
poplist <- popTableAsList(popids)
if(verbose) print(paste("Read population IDs for",length(poplist),"populations"))


popdata <- t( matColMeans(t(data),poplist))

write.table(popdata, file=files[3], quote=FALSE, row.names=TRUE, col.names=TRUE)
if(verbose) print(paste("Wrote reference panel to",files[3]))
