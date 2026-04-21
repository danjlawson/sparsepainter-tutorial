## Runs admixture estimation on the reference panels

if(!exists("args") || class(args)!="character"){
    args <- commandArgs(trailingOnly = TRUE)
}

if(length(args)<1){
    cat("Usage: Rscript admixture_vs_panel.R <options> <pop ids> <admix results> <output root>\n")
    cat("OPTIONS:\n")
    cat("-v : Verbose mode (default: FALSE)\n")
    cat("Example usage: Rscript ../code/panel_confusion_matrix.R -v testcombined.pop.ids results/test_admix_vs_panel.txt results/confusion")
    stop("Must provide 1 or more arguments")
}

i=1
tempargs=args
verbose=FALSE

while(i<=length(tempargs)) {
  if(tempargs[i]=="-v"){
        verbose=TRUE
        tempargs=tempargs[-i]
    } else {
        i=i+1
    }
}
files = tempargs
if(length(files)!=3){
    stop("Must provide exactly 3 file names: population_ids, admix_results, and output root")
}

source("../code/FinestructureLibrary.R") # read in the R functions, which also calls the needed packages

idtable = read.table(files[1], header=TRUE, row.names=1)
if(verbose) print(paste("Read population IDs for",nrow(idtable),"individuals and",length(unique(idtable[,1])),"populations"))
idlist = popTableAsList(idtable)

admix <- as.matrix(read.table(files[2], header=TRUE, row.names=1))
if(verbose) print(paste("Read admixture panel with",nrow(admix)," surrogate populations and",ncol(admix),"donor populations"))

confusion = t(matColMeans(t(admix),idlist))
if(all(rownames(confusion)%in%colnames(confusion))){
    confusion = confusion[,rownames(confusion)]
}

write.table(confusion, file=paste0(files[3],".confusion.csv"), quote=FALSE, row.names=TRUE, col.names=TRUE)
if(verbose) print(paste("Wrote confusion table to",paste0(files[3],".confusion.csv")))

require("pheatmap")
## Make a heatmap of the confusion matrix

## These parameters might need to be adjusted for larger or smaller datasets, but should work well for a few dozen populations
## Increase the "number_format" parameter to show more decimal places in the confusion matrix values, useful for smaller datasets
pdf(paste0(files[3],".confusion.pdf"), width=20, height=20)
pheatmap(confusion, cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=TRUE, number_format="%.1f", main="Admixture panel confusion matrix")
dev.off()
if(verbose) print(paste("Wrote confusion heatmap to",paste0(files[3],".confusion.pdf")))