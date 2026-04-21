## Runs admixture estimation on the reference panels

if(!exists("args") || class(args)!="character"){
    args <- commandArgs(trailingOnly = TRUE)
}

if(length(args)<1){
    cat("Usage: Rscript admixture_vs_panel.R <options> <admix panel> <target chunklengths> <output file>\n")
    cat("OPTIONS:\n")
    cat("-v : Verbose mode (default: FALSE)\n")
    cat("Example usage: Rscript ../code/admixture_vs_panel.R -v test_admix_panel.txt sparsepaint/test.combined.chunklength.txt results/test_admix_vs_panel.txt\n")
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
    stop("Must provide exactly 3 file names: admix_panel, targetfile, and output")
}

source("../code/FinestructureLibrary.R") # read in the R functions, which also calls the needed packages

admix_panel <- as.matrix(read.table(files[1], header=TRUE, row.names=1))
if(verbose) print(paste("Read admixture panel with",nrow(admix_panel)," surrogate populations and",ncol(admix_panel),"donor populations"))
targetfile <- as.matrix(read.table(files[2], header=TRUE, row.names=1))
if(verbose) print(paste("Read target file with",nrow(targetfile),"individuals and",ncol(targetfile),"donor populations"))

require(nnls)

## General NNLS functions

admix.nnls<-function(X,Y){
    ## Do an nnls admixture analysis (X ~ Y) where X is a vector
    if(!"numeric"%in% class(X)) stop("X must be numeric in admix.nnls")
    ourmix=getoverallfit(Y,X)$x
    ourmix=ourmix/sum(ourmix)
    ourmix
}

admix.nnls.all<-function(X,Y,verbose=TRUE){
    ## Do nnls admixture on each ROW of the MATRIX x
    if(!"matrix" %in% class(X)) stop("X must be a matrix in admix.nnls.all")
    
    ret<-t(sapply(1:dim(X)[1],function(i){
        if(verbose) print(paste("Processing ind",i,"of",dim(X)[1]))
        admix.nnls(X[i,],Y)
    }))
    rownames(ret)<-rownames(X)
    ret
}


#######################
## MAIN NNLS FUNCTIONS, from GLOBETROTTER paper

getfit=function(predmat,fitdata,restrict=1){
    ## Gets the nnls fit with the restriction that a given row is not used.
    ## This addresses the sum to one constraint.
    ## If we succeed in getting a valid return, it is also guarenteed to be the
    ## best fitting. Also, we are guaranteed to get this from some restriction
    ## choice
  temp=predmat[-restrict,,drop=FALSE]
  for(i in 1:nrow(temp)) temp[i,]=temp[i,]-predmat[restrict,]

  fitdata2=fitdata-predmat[restrict,]

  v=nnls(t(temp),fitdata2)
  
  x=v$x
  newx=1:nrow(predmat)
  newx[!((1:nrow(predmat))==restrict)]=x
  newx[restrict]=1-sum(x)
  v$x=newx
  names(v$x)=rownames(predmat)

  return(v)
}

getoverallfit=function(predmat,fitdata){
  restrict=1
  rep=1
  i=1
  while(rep==1){
    q=getfit(predmat,fitdata,restrict=i)

    if(q$x[i]>0) rep=0
    i=i+1
  }

  return(q)
}

admix=admix.nnls.all(targetfile,admix_panel,verbose=verbose)
write.table(admix, file=files[3], quote=FALSE, row.names=TRUE, col.names=TRUE)
if(verbose) print(paste("Wrote admixture results to",files[3]))