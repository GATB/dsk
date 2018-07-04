#!/usr/bin/Rscript

# Author : Claire Lemaitre


#usage :  ./plot-histo2D.R  dsk.histo2D
# or      RScript plot-histo2D.R  dsk.histo2D
args=commandArgs(TRUE)
if (length(args)<1) {
  cat("Usage:\n")
  cat("./plot-histo2D.R dsk.histo2D [xmax]\n")
  cat("  1 obligatory argument :\n      - kmer histo2D file (as output by DSK with option -histo2D)\n")
  cat("  optional arguments :\n      - max x value to show on the x axis (default =  automatically inferred)\n")
  cat("  output :  dsk.histo2D.png (png image file)\n")
  cat("  note :  plot rendering is inspired from KAT (y axis is not in log scale)\n")
  quit()
}


hist.file=args[1]

tab=read.table(hist.file)


mat = tab[,-1]

# automatic adjustement of ymax :
#first compute the kmer profile of the read dataset (sum over all columns, ie. all assembly abundances)
linetot=apply(mat[-c(1,nrow(mat)),],1,'sum')
#eliminating the last line of mat because cumul of all abundances > abundance_max (def = 10001)
# then identifying the first increase :
beg=which(diff(linetot)>0)[1]
# then finding the maximum after the first increase (hopefully skipping the highly numerous kmers with low abundance and containing sequencing errors)
ymax=max(linetot[beg:length(linetot)])*1.05

#ymax=max(mat[,2])*1.05

if (length(args)>1) {
  xmax=as.numeric(args[2])
}else{
  #auto computed as the largest abundance with kmer count > 0.5*ymax/100  (0.5%)
  xmax=max(which(linetot>=0.5*ymax/100))
  #print(xmax)
}

if(ymax>1e+6){
  mat=mat/1e+6
  ymax=ymax/1e+6
  unit="1e+6"
}else{
  mat=mat/1e+3
  ymax=ymax/1e+3
  unit="1e+3"
}

mycolors=c("black","red","mediumpurple2","palegreen2","steelblue3","peachpuff1")

png.file=paste(hist.file,".png",sep="")
bitmap(png.file,"png256",width=7,height=6,res=300)
suppressWarnings (barplot(t(mat[1:xmax,c(1:6)]),beside=F,ylim=c(0,ymax),col=mycolors,border=NA,space=0,axes=F,axisnames=F,ylab=paste("Number of distinct kmers (",unit,")",sep=""),xlab="kmer multiplicity",xpd=F,main="kmer comparison plot"))
axis(1,at=pretty(c(0,xmax)))#+1,labels = pretty(c(0,xmax)) )
axis(2,at=pretty(c(0,ymax)))
grid(lty=1)
legend(xmax,ymax,fill=mycolors,legend=paste(0:5,"x",sep=""),xjust=1)
d=dev.off()

cat(paste("... done, image output in file",png.file,"\n"))

