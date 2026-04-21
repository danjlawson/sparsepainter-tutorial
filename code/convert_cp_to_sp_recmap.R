## Converts recombination maps from chromopainter to sparsepainter format

if(!exists("args") || class(args)!="character"){
    args <- commandArgs(trailingOnly = TRUE)
}

if(length(args)<1){
    cat("Usage: Rscript estimate_lambda.R <options> <files>\n")
    cat("OPTIONS:\n")
    cat("-i <file>: Input file ending (default: .recombfile)\n")
    cat("-o <file>: Output file ending, which will replace the input file ending (default: .map)\n")
    cat("-v : Verbose mode (default: FALSE)\n")
    cat("Example usage: Rscript ../code/convert_cp_to_sp_recmap.R -v rawdata/EuropeSample.small.chrom{1..22}.recombfile
\n\n")
    stop("Must provide 1 or more arguments")
}

i=1
inending=".recombfile"
outending=".sp.map"
tempargs=args
verbose=FALSE
while(i<=length(tempargs)) {
    if(tempargs[i]=="-o"){
        i=i+1
        outending=tempargs[i]
        tempargs=tempargs[-c(i-1,i)]
    } else if(tempargs[i]=="-i"){
        i=i+1
        inending=tempargs[i]
        tempargs=tempargs[-c(i-1,i)]
    } else if(tempargs[i]=="-v"){
        verbose=TRUE
        tempargs=tempargs[-i]
    } else {
        i=i+1
    }
}
files=tempargs

for(f in files){
    outf=sub(inending, outending, f)
    if(outf==f){
        stop("ERROR: input file ending not found in file name. Output file would have the same name as input file.\n")
    }
    if(verbose){
        cat("Processing file:", f, "->", outf, "\n")
    }
    rin=read.table(f,header=T)
    rin$diffbp=c(0,diff(rin[,1]))
    rin$recom.diff=c(0,rin[-1,"diffbp"]*rin[-dim(rin)[1],"recom.rate.perbp"]*100) # 100: convert from Morgans to centi-Morgans
    rin$cum.recom.cm=cumsum(rin$recom.diff)
    rout=rin[,c("start.pos","cum.recom.cm")]
    colnames(rout)<-c("position(base)", "genetic_distance(centiMorgan)")
    write.table(rout,file=outf,row.names=FALSE,quote=FALSE)
}
