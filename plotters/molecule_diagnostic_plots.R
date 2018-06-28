#! /usr/bin/env Rscript

## molecule_diagnostics.R
## Eva KF Chan
## https://github.com/ekfchan
## Script generates several diagnostic plots related to the molecule length and label density, partitioned by mapping "batches".
## The input file is derived from BNX v1.3 (e.g. all_sorted.bnx) by pulling out all 0-channel entries: grep "^0" all_sorted.bnx > in.txt 

## Arguments: 
## ==========
## Requires: 
## 	channel0:  channel 0 entries from bnx v1.3 file; e.g. dat.bnx
##	out.pdf: Filename for the PDF output of this script. 

args <- commandArgs(TRUE)
if( length( args )!=2 ) { stop("usage: ./molecule_diagnostics.R in.txt out.pdf\n") }

infile = args[1]
outfile = args[2]
cat( "BNX:", infile, "\n" )
cat( "OUTPUT:", outfile, "\n" )

if(!file.exists(infile)) { stop("Can't file",infile,"\n") }
if(file.exists(outfile)) { stop(outfile,"already exists!\n") }

## Needed Functions: 
## =================
wrap.labels <- function(x, len)
{
  sapply(lapply(strwrap(gsub(""," ",x),len*2,simplify=F),gsub,pattern="[[:space:]]",replacement=""),paste,collapse="\n")
}

plotboxes <- function( fac, facname ) { 
	par(mfrow=c(2,2), oma=c(5,5,3,5), mar=rep(0,4)) 
	# Molecule Length
	znames = unique(fac)
	boxplot( log10(dat$Length) ~ fac, las=1, axes=F, ann=F, las=1, outline=F )
	axis( 2, at=log10(c(100000,150000,300000,500000,1000000,1500000,2000000)), labels=c("100 kb","150 kb",paste(c(0.3,0.5,1.0,1.5,2.0),"Mb",sep="")), lwd=0, lwd.ticks=1, las=1 )
	mtext( text="Molecule Length (Mb)", side=2, line=3 )
	box()
	# Label Density 
	boxplot( (dat$NumberofLabels/dat$Length*100000) ~ fac, ann=F, las=1, axes=F, outline=F )
	axis( 4, las=1, lwd=0, lwd.ticks=1 )
	mtext( text="Label Density (/ 100 Mb)", side=4, line=3 )
	box()
	# Average Intensity 
	boxplot( log10(dat$AvgIntensity) ~ fac, ann=F, axes=F, las=1, outline=F )
	axis( 2, lwd=0, lwd.ticks=1, las=1, at=c(2,log10(500),3,log10(5000),4,log10(50000),5), labels=c(100,500,expression(10^3),"",expression(10^4),"",expression(10^5)) )
	mtext( text="Average Intensity", side=2, line=3 )
	axis( 1, lwd=0, lwd.ticks=1, las=1, at=1:length(znames), labels=wrap.labels(znames,16), cex.axis=0.8, padj=1, line=0 )
	box()
	# SNR 
	boxplot( log10(dat$SNR) ~ fac, ann=F, axes=F, las=1, outline=F )
	axis( 4, lwd=0, lwd.ticks=1, las=1, at=0:4, labels=c(1,10,100,expression(10^3),expression(10^4)) )
	mtext( text="SNR", side=4, line=3 )
	axis( 1, lwd=0, lwd.ticks=1, las=1, at=1:length(znames), labels=wrap.labels(znames,16), cex.axis=0.8, padj=1, line=0 )
	box()
	# Title
	mtext( text=facname, side=3, line=1, outer=T )
}

plotlines <- function( fac, facname ) { 
	par(mfrow=c(2,2), oma=c(5,5,3,5), mar=rep(0,4)) 
	# Molecule Length
	z <- split(log10(dat$Length),f=fac)
	plot( x=as.numeric(names(z)), y=sapply(z,mean), las=1, ann=F, axes=F, las=1, type="l", ylim=log10(c(100000,1000000)), lwd=2)
	lines( as.numeric(names(z)), sapply(z,quantile,0.25), lty=1, col="grey" )
	lines( as.numeric(names(z)), sapply(z,quantile,0.75), lty=1, col="grey" )
	axis( 2, at=log10(c(100000,150000,300000,500000,1000000)), labels=paste(c(0.10,0.15,0.30,0.5,1.0),"Mb",sep=""), lwd=0, lwd.ticks=1, las=1 )
	mtext( text="Molecule Length (Mb)", side=2, line=3 )
	box()
	rm(z)
	# Label Density 
	z <- split( (dat$NumberofLabels/dat$Length*100000), f=fac)
	plot( x=as.numeric(names(z)), y=sapply(z,mean), las=1, ann=F, axes=F, las=1, type="l", ylim=c(5,25), lwd=2 )
	lines( as.numeric(names(z)), sapply(z,quantile,0.25), lty=1, col="grey" )
	lines( as.numeric(names(z)), sapply(z,quantile,0.75), lty=1, col="grey" )
	axis( 4, las=1, lwd=0, lwd.ticks=1 )
	mtext( text="Label Density (/ 100 Mb)", side=4, line=3 )
	box()
	rm(z)
	# Average Intensity 
	z <- split( log10(dat$AvgIntensity), f=fac)
	plot( x=as.numeric(names(z)), y=sapply(z,mean), las=1, ann=F, axes=F, las=1, type="l", ylim=c(2,5), lwd=2 )
	lines( as.numeric(names(z)), sapply(z,quantile,0.25), lty=1, col="grey" )
	lines( as.numeric(names(z)), sapply(z,quantile,0.75), lty=1, col="grey" )
	axis( 2, lwd=0, lwd.ticks=1, las=1, at=2:5, labels=c(100,expression(10^3),expression(10^4),expression(10^5)) )
	mtext( text="Average Intensity", side=2, line=3 )
	axis( 1, lwd=0, lwd.ticks=1, las=1 )
	box()
	rm(z)
	# SNR 
	z <- split( log10(dat$SNR), f=fac)
	plot( x=as.numeric(names(z)), y=sapply(z,mean), las=1, ann=F, axes=F, las=1, type="l", ylim=c(0,4), lwd=2 )
	lines( as.numeric(names(z)), sapply(z,quantile,0.25), lty=1, col="grey" )
	lines( as.numeric(names(z)), sapply(z,quantile,0.75), lty=1, col="grey" )
	axis( 4, lwd=0, lwd.ticks=1, las=1, at=0:4, labels=c(1,10,100,expression(10^3),expression(10^4)) )
	mtext( text="SNR", side=4, line=3 )
	axis( 1, lwd=0, lwd.ticks=1, las=1 )
	box()
	rm(z)
	# Title
	mtext( text=facname, side=3, line=1, outer=T )
}


## Data: 
## =====
# Test file for header
mycon <- file(infile, open="r")
dathead <- readLines(mycon, n=1, warn=F)
hasHead=F
if( length(grep("^0",dathead))==0 | length(grep("^#",dathead))>0 ) { hasHead=T }
if( hasHead ) { dathead = strsplit(dathead[1],split="\t")[[1]] } else { dathead=c("LabelChannel","MoleculeID","Length","AvgIntensity","SNR","NumberofLabels","OriginalMoleculeId","ScanNumber","ScanDirection","ChipId","Flowcell","RunId","Column","StartFOV","StartX","StartY","EndFOV","EndX","EndY","GlobalScanNumber") } 

# Load data 
dat <- read.table( file=mycon, header=hasHead, sep="\t", stringsAsFactors=F, na.strings="", col.names=dathead )
close(mycon)

dat$LabelChannel <- NULL


## Make Diagnostic Plots: 
## ======================
pdf(file=outfile, width=10, height=10)

## ~~ Overview ~~  
par(mfrow=c(2,2)) 
# Molecule length density
z <- density( log10(dat$Length) )
plot(x=log10(c(100000,2000000)), y=c(0,ceiling(max(z$y))), main="Filtered Molecule Length", las=1, axes=F, pch=NA, sub=paste("N =",format(z$n,big.mark=",")," Bandwidth =",round(z$bw,6)), ylab="Density", xlab="" )
abline( v=log10(c(100000,150000,300000,500000,1000000,1500000,2000000)), col="grey" )
lines( z, lwd=2 )
axis( 1, at=log10(c(100000,150000,300000,500000,1000000,1500000,2000000)), labels=paste(c(0.10,0.15,0.3,0.5,1.0,1.5,2.0),"Mb"), lwd=0, lwd.ticks=1, cex.axis=0.8 )
axis(2,las=1,lwd.ticks=1,lwd=0)
box()
# Label density 
plot(density(dat$NumberofLabels/dat$Length*100000), las=1, main="Label Density (/ 100 kb)")
abline(v=mean(dat$NumberofLabels/dat$Length*100000), lty=3 )
# Intensity
plot(density(log10(dat$AvgIntensity)), las=1, main="Average Intensity", axes=F, xlim=c(2,5))
axis(2,las=1,lwd.ticks=1,lwd=0)
axis(1,las=1,at=2:5,labels=c(100,"1,000",expression(10^4),expression(10^5)),lwd=0,lwd.ticks=1)
box()
# SNR
plot(density(log10(dat$SNR)), las=1, main="SNR", axes=F, xlim=c(0,4))
axis(2,las=1,lwd.ticks=1,lwd=0)
axis(1,las=1,at=0:4,labels=c(0,10,100,expression(10^3),expression(10^4)),lwd=0,lwd.ticks=1)
box()

## ~~ Split by Factors ~~ 
for( facname in c("ScanNumber","ScanDirection","ChipId","Flowcell","RunId","Column","StartFOV","GlobalScanNumber") ) {
	fac = dat[,facname]
	if( length(unique(fac))>1 ) {
		if( length(unique(fac))<10 ) {
			plotboxes(fac=fac, facname=gsub("([a-z])([A-Z])","\\1 \\2",facname))
		} else {
			plotlines(fac=fac, facname=gsub("([a-z])([A-Z])","\\1 \\2",facname))
		}
	}
	rm( fac )
} 

dev.off()

q('no')


