#! /usr/bin/env Rscript

## wgCNV.R
## Eva KF Chan
## https://github.com/ekfchan
## Script to plot genome-wide CNV profile, taken as the fractionalCopyNumber from output/contigs/alignmolvref/copynumber/cnv_rcmap_exp.txt of the Bionano Solve v3 pipeline output

## Arguments: 
## ==========
## Inputs: 
##	- output/contigs/alignmolvref/copynumber/cnv_rcmap_exp.txt
## Output: 
##	- PDF

args <- commandArgs(TRUE)
if( length( args )!=2 ) { stop("usage: ./wgCNV.R cnv_rcmap_exp.txt out.pdf\n") }

infile = args[1]
if( !file.exists(infile) ) { stop("Cannot find ",infile,"\n") }
outfile = args[2]


## Setup:
## ======
# Define UCSC default chromosome colours [original colour definition taken from CIRCOS (https://github.com/vigsterkr/circos/blob/master/tutorials/6/10/colors.conf) as RGB and converted to hexadecimal string]
chrcolors <- strsplit("#996600,#666600,#99991E,#CC0000,#FF0000,#FF00CC,#FFCCCC,#FF9900,#FFCC00,#FFFF00,#CCFF00,#00FF00,#358000,#0000CC,#6699FF,#99CCFF,#00FFFF,#CCFFFF,#9900CC,#CC33FF,#CC99FF,#666666,#999999,#999999,#CCCCCC,#CCCCCC,#CCCC99,#CCCC99,#79CC3D,#FFFFFF",",")[[1]]
# chrcolors = paste(chrcolors,"5A",sep="")	#add 90/255 transparency 
names(chrcolors) <- as.character(c(1:23,"X","24","Y","M",0,"Un",NA))

## Data: 
## =====
rcmap <- read.table(file=infile, header=T, comment.char="", sep="\t", stringsAsFactors=F, na.strings="" )
names(rcmap) = sub("X\\.","",names(rcmap))	#header line generally starts with "#" symbol which will be imported into R as "X."

yrange = range(rcmap$fractionalCopyNumber)
yrange = c(floor(yrange[1]),ceiling(yrange[2]))
xmax=sum(rcmap$CMapId>=1 & rcmap$CMapId<=22)+1
gridelines = c(yrange[1],0,2,seq(5,yrange[2],by=5))

pdf( file=outfile, width=11, height=6 ) 
par(mar=c(4,4,2,2), mgp=c(2,1,0))

plot( x=c(0,xmax), y=yrange, pch=NA, axes=F, ann=F )
segments(x0=rep(0,length(gridelines)),x1=rep(xmax,length(gridelines)), y0=gridelines,y1=gridelines,col="grey")
for( i in 1:22 ) {	#Bionano advices that CNV estimates for sex chromosomes are not reliable.
	inds <- which(rcmap$CMapId==i)
	lines( x=inds, y=rcmap$fractionalCopyNumber[inds], col=chrcolors[i], lwd=2 )
	rm(inds)
}

xmarks <- c(match(1:22,rcmap$CMapId),xmax)
axis( 1, at=xmarks, labels=F, lwd=1, lwd.ticks=1, pos=yrange[1]-1 )
axis( 1, at=xmarks[-length(xmarks)]+((xmarks[-1]-xmarks[-length(xmarks)])/2), labels=1:22,tick=F,lty=0, cex=0.9)
mtext( text="Chr:", side=1, line=1, at=yrange[1], adj=1 )
axis( 2, at=gridelines, las=1, pos=0, lwd=1, lwd.ticks=1 )
mtext( text="Copy Number", side=2, line=1 )

# dev.off()
graphics.off()


q('no')
