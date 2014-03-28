#!/usr/bin/Rscript


pdf("figure.pdf",width=5,height=3)
data <-read.csv(file="ecoli_log",sep="\t",head=TRUE)


layout(matrix(c(1,2,3,3), 2, 2, byrow = T), 
       widths=c(0.4,0.4), heights=c(0.8,0.15))


par(mar=c(4,4,3,2.3)) # adjust margins
# transform times to seconds
parse_time <- function(s) as.integer(strsplit(s, "m")[[1]][1])*60+as.integer(strsplit(strsplit(s, "m")[[1]][2],"[.]")[[1]][1])
parse_times <- function(s) sapply(as.character(s),parse_time)
times <- parse_times(data$Times)

# plot x axis = memory, y axis = time, 1 line per disk usage
disk <- unique(data$Disk)
disk_tics <- c(1,2,3,4,5,6)
memory_values <- c(10,100,1000)
plot(range(disk_tics),range(c(0,500)),main="E. coli DNA",xlab="Disk space (MB)",ylab="Time (s)",xaxt="n",type="n")
axis(1, at=1:6, labels=disk)
j <- 1
#for (i in unique(data$Memory)) {  
for (i in memory_values) {  
  lines(disk_tics,parse_times(subset(data,Memory==i)$Time),col=4-j,lwd=2,lty=4-j)
  j <- j+1
}

par(mar=c(4,4,3,2.3)) # adjust margins

# droso plot
data_droso <- read.csv(file="droso_log",sep="\t",head=TRUE)
times <- parse_times(data_droso$Times)
disk <- unique(data_droso$Disk)
plot(range(disk_tics),range(c(0,1500)),main="Drosophila RNA",xlab="Disk space (MB)",ylab="Time (s)",xaxt="n",type="n")
axis(1, at=1:6, labels=disk)
j <- 1
for (i in memory_values) {  
  lines(disk_tics,parse_times(subset(data_droso,Memory==i)$Time),col=4-j,lwd=2,lty=4-j)
  j <- j+1
}

# legend common to both plots
par(mar=c(0,0,0,0)) # adjust margins and remove borders
plot(1, type = "n", axes=FALSE, xlab="", ylab="", bty="n", frame.plot=FALSE) # just a transparent plot
legend("top", legend=memory_values, cex=1,lty=c(3,2,1),col=c(3,2,1),title="Memory (MB)", lwd=2, horiz = TRUE, box.col = "white" );

