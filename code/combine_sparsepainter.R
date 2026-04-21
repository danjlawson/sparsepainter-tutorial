## Converts recombination maps from chromopainter to sparsepainter format

if(!exists("args") || class(args)!="character"){
    args <- commandArgs(trailingOnly = TRUE)
}

if(length(args)<1){
    cat("Usage: Rscript combine_sparsepainter.R <options> <files>\n")
    cat("OPTIONS:\n")
    cat("-o <file>: Output file ending\n")
    cat("-v : Verbose mode (default: FALSE)\n")
    cat("-c <n>: Number of times to expect each individual in the input file (default: 2). If >1, looks for form _<i-1> at the end of individual IDs to identify replicates. SparsePainter produces -n 2 format data for diploids.\n\n")
    cat("Example usage: Rscript ../code/convert_cp_to_sp_recmap.R -v -o sparsepaint/test.combined.chunklength.txt sparsepaint/test{1..22}.sp.loo_chunklength.txt.gz
\n\n")
    stop("Must provide 1 or more arguments")
}

i=1
outfile="combine_sparsepainter_results.txt"
tempargs=args
verbose=FALSE
copies=2

while(i<=length(tempargs)) {
    if(tempargs[i]=="-o"){
        i=i+1
        outfile=tempargs[i]
        tempargs=tempargs[-c(i-1,i)]
    } else if (tempargs[i]=="-c"){
        i=i+1
        copies=as.numeric(tempargs[i])
        tempargs=tempargs[-c(i-1,i)]
    } else if(tempargs[i]=="-v"){
        verbose=TRUE
        tempargs=tempargs[-i]
    } else {
        i=i+1
    }
}

files=tempargs

for (f in files){
    if(verbose) print(paste("Reading",f))
    rin=read.table(f,header=T,row.names=1)
    if(f == files[1]){
        res = rin
    } else {
        res = res + rin
    }
}

data = res
if (copies>1){
    if(verbose) print(paste("Looking for",copies,"copies of each individual in the input file"))
    datalist <- lapply(1:copies, function(i) {
        if(verbose) print(paste("Processing copy",i))
        copydata <- data[grep(paste0("_",i-1,"$"), rownames(data)), , drop=FALSE]
        rownames(copydata) <- sub(paste0("(_",i-1,")$"), "", rownames(copydata))
        copydata
    })
    data=Reduce("+", datalist)
}

write.table(data,file=outfile,row.names=TRUE,quote=FALSE)
if(verbose) print(paste("Written combined results to",outfile))