#!/usr/bin/Rscript

# Author : Claire Lemaitre


#usage :  ./plot-histo.R  dsk.histo
# or      RScript plot-histo.R  dsk.histo
args=commandArgs(TRUE)
if (length(args)<1) {
cat("Usage:\n")
cat("./plot-histo.R dsk.histo [xmax]\n")
cat("  1 obligatory argument :\n      - kmer histo file (as output by DSK)\n")
cat("  optional arguments :\n      - max x value to show on the x axis (default =  max in the histo file)\n")
cat("  output :  dsk.histo.png (png image file)\n")
cat("  note :  the y axis is in log scale\n")

quit()
}


hist.file=args[1]


tab=read.table(hist.file)[,1:2]

if (length(args)>1) {
    xmax=as.numeric(args[2])
}else{
    xmax=max(tab$V1)
}

png.file=paste(hist.file,".png",sep="")
bitmap(png.file,"png256",width=7,height=6,res=300)
suppressWarnings (plot(tab$V1,tab$V2,xlim=c(0,xmax),type="n",log="y",xlab="Kmer abundance",ylab="Number of distinct kmers",main="Kmer profile"))
grid(lty=1)
lines(tab$V1,tab$V2)
d=dev.off()

cat(paste("... done, image output in file",png.file,"\n"))
