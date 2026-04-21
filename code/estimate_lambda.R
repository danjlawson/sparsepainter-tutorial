

## Description: Script to combine lambda estimates from separate SparsePainter runs and create a file containing their average for SparsePainter use.
## Author: Daniel J Lawson
## Date: Feb 2026
## Email: dan.lawson@bristol.ac.uk
## License: GPL V3

if(!exists("args") || class(args)!="character"){
    args <- commandArgs(trailingOnly = TRUE)
}

if(length(args)<1){
    cat("Usage: Rscript estimate_lambda.R <options> <files>\n")
    cat("OPTIONS:\n")
    cat("-o <file>: Output file for the estimated lambda (default: lambda.txt)\n")
    cat("-v : Verbose mode (default: FALSE)\n\n")
    cat("Example usage: Rscript ../code/estimate_lambda.R -v -o sparsepaint/lambda.txt sparsepaint/test{1..22}.sp.estlambda_fixedlambda.txt\n\n")
    stop("Must provide 1 or more arguments")
}

i=1
outlambda="lambda.txt"
verbose=FALSE
tempargs=args
while(i<=length(tempargs)) {
    if(tempargs[i]=="-o"){
        i=i+1
        outlambda=tempargs[i]
        tempargs=tempargs[-c(i-1,i)]
    } else if(tempargs[i]=="-v"){
        verbose=TRUE
        tempargs=tempargs[-i]
    } else {
        i=i+1
    }
}

files=tempargs
## Read in each file in files, extract the lambda estimate, and average them
## They are separated by a : and the lambda estimate is the second element
lambdas=c()
for(file in files){
    dat=read.table(file,header=FALSE,stringsAsFactors=FALSE, sep=":")
    print(dat)
    lambdas=c(lambdas, as.numeric(dat[2]))
    if(verbose){
        cat(paste("Read lambda estimate", dat[2], "from file", file, "\n"))
    }
}
lambda=mean(lambdas)
if(verbose){
    cat("Average lambda estimate across files:", lambda, "\n")
    cat("Writing average lambda estimate to file:", outlambda, "\n")
}
write.table(lambda,file=outlambda,row.names=FALSE,col.names=FALSE,quote=FALSE)