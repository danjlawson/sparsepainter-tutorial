## --- Data ---
## input data
for(chr in 1:22){
  
datatext <-readLines(paste0("targetpaint/test",chr,".sp.local_prob.txt.gz"))
## Extract the column headers, start and end row IDs for each haplotype
colnames=strsplit(datatext[2]," ")[[1]]
colnames[1:2]<-c("start","end")
splitpts=which(sapply(strsplit(datatext," "),length)==1)
endpts=c(splitpts[-1]-1,length(datatext))
hapnames=as.character(strsplit(datatext[splitpts]," "))

## Extract the data as a matrix
splitdf=as.data.frame(t(data.frame(start=splitpts,end=endpts)))
datalist=lapply(splitdf,function(x){
    tmp=do.call("rbind",strsplit(datatext[x[1]:x[2]]," ")[-1])
    tmp=apply(tmp,2,as.numeric)
    colnames(tmp)=colnames
    tmp
})
## Ensure that no rounding errors lead to sums that are less than one
datalist = lapply(datalist,function(data){
    data[,-(1:2)] <- data[,-(1:2)]/rowSums(data[,-(1:2)])
    as.data.frame(data)
})

# --- Colours for each ancestry class ---
ancestry_cols <- c("#E41A1C", "#377EB8", "#4DAF4A")  # red, blue, green


png(paste0("results/localAncestry_chr",chr,".png"),height=1500,width=1000)
layout(matrix(1:(length(datalist)+1),ncol=1),
       heights=c(rep(1,length(datalist)),1)
       )
par(mar=c(1,4,1,4),cex=2)
for(i in 1:length(datalist)){
    data=datalist[[i]]
    if(i==1) main="Local Ancestry"
    else main=""
    if(i==length(datalist)) xlab="SNP number"
    else xlab=""
    ## --- Plot ---
    ## Use the full genomic range as x-axis limits
    x_min <- min(data$start)
    x_max <- max(data$end)

    ## Empty plot canvas
    plot(NULL,
         xlim = c(x_min, x_max),
         ylim = c(0, 1),
         xlab = xlab,
         ylab = "Ancestry proportion",
         main = main,
     xaxs = "i", yaxs = "i")  # no padding at axis edges

    ## Draw each segment as stacked filled rectangles
    for (i in seq_len(nrow(data))) {
        seg_start <- data$start[i]
        seg_end   <- data$end[i] + 1  # extend to the start of the next position
        probs     <- as.numeric(data[i, -(1:2)])
        
        ## Compute cumulative heights for stacking
        cum_probs <- c(0, cumsum(probs))
        
        for (k in seq_along(probs)) {
            if (probs[k] > 0) {
                rect(
                    xleft   = seg_start,
                    xright  = seg_end,
                    ybottom = cum_probs[k],
                    ytop    = cum_probs[k + 1],
                    col     = ancestry_cols[k],
                    border  = NA  # no borders → smooth appearance
                )
            }
        }
    }
}
  plot(NULL,
         xlim = c(x_min, x_max),
         ylim = c(0, 1),
         xlab = "",
         ylab = "",
         main = "",
       xaxs = "i", yaxs = "i",axes=FALSE)  # no padding at axis edges
legend("center",legend=colnames[-(1:2)],fill=ancestry_cols,cex=0.75,bty="n")
dev.off()

}
