#! /usr/bin/env Rscript

## bionano_ideogram.R
## Eva KF Chan
## https://github.com/ekfchan
## Script to plot the alignment coverage of a genome map set in the cytoband ideogram format.

## Arguments: 
## ==========
## Requires: 
## 	ideogram.txt:  This is a Cytoband Ideogram file, which can be obtained from the collection of UCSC Tables. e.g. cytoBandIdeo table from hg38 database last updated 2014-06-11 (see data/hg38_goldenpath_cytoBandIdeo.txt)
##	bionano_xmap.bed:  This is essentially columns 3 (RefContigID), 6 (RefStartPos), and 7 (RefEndPos) of one of the xmap outputs from the Bionano Solve pipeline: e.g. output/contigs/exp_refineFinal1_sv/EXP_REFINEFINAL1_full.xmap
##	out.pdf: Filename for the PDF output of this script. 

args <- commandArgs(TRUE)
if( length( args )!=3 ) { stop("usage: ./bionano_ideogram.R ideogram.txt bionano_xmap.bed out.pdf\n") }

ideofile = args[1]
xmapfile = args[2]
outfile = args[3]
cat( "Ideogram:", ideofile, "\n" )
cat( "XMAP:", xmapfile, "\n" )
cat( "OUTPUT:", outfile, "\n" )

## Ideogram setup
## ==============
# "Darker bands are AT-rich, and lighter bands are GC-rich. Moving the mouse pointer over the cytobands will show the cytoband type. 'gneg' and 'gpos' define the stain intensity, 'acen' means acrocentric bands, 'gvar' stands for variable heterochromatic region and 'stalk' refers to tightly constricted regions on the short arms of acrocentric chromosomes." [https://www.genomatix.de/online_help/help_eldorado/GenomeBrowser.html]
# karyotype colors were based on those used by CIRCOS [https://github.com/vigsterkr/circos/blob/master/tutorials/6/10/colors.conf], with modifications made to gvar = pink & stalk = powderblue
gieColor = scan(what="character", sep="\n", quiet=T)
gpos100 = 0,0,0
gpos    = 0,0,0
gpos75  = 130,130,130
gpos66  = 160,160,160
gpos50  = 200,200,200
gpos33  = 210,210,210
gpos25  = 200,200,200
# gvar    = 220,220,220
gvar    = 255,192,203
gneg    = 255,255,255
acen    = 217,47,39
# stalk   = 100,127,164
stalk   = 176,224,230



gieColor=strsplit(gieColor,split="\\s+=\\s+")
gieColor <- data.frame(name=unlist(lapply(gieColor,'[',1)), rgb=unlist(lapply(gieColor,'[',2)), stringsAsFactors=F)
gieColor$red = as.integer(sub("^(.+),(.+),(.+)$","\\1",gieColor$rgb))
gieColor$green = as.integer(sub("^(.+),(.+),(.+)$","\\2",gieColor$rgb))
gieColor$blue = as.integer(sub("^(.+),(.+),(.+)$","\\3",gieColor$rgb))
gieColor$hex = rgb(red=gieColor$red, green=gieColor$green, blue=gieColor$blue, names=gieColor$name, maxColorValue=255)
 
## Ideogram Data:
## ==============
cat( "Loading ideogram:", ideofile, "\n")
ideo <- read.table( file=ideofile, header=F, stringsAsFactors=F, na.strings="", sep="\t", col.names=c("chrom","chromStart","chromEnd","name","gieStain") )
ideo <- ideo[-grep("_",ideo$chrom),]	#exclude alt and random chrs
ideo <- ideo[which(ideo$chrom!="chrM"),]	#exclude chrM
ideo$chr = as.integer(sub("chr","",sub("Y",24,sub("X",23,ideo$chrom))))


## Xmap Data: 
## ==========
# xmap <- read.table( file=xmapfile, header=F, sep="\t", stringsAsFactors=F, col.names=c("chrom","start","end","counts") )
# xmap$height <- xmap$counts; xmap$height[which(xmap$counts>100)]=100
xmap <- read.table( file=xmapfile, header=F, sep="\t", stringsAsFactors=F, col.names=c("chrom","start","end") )

## Plot
## ====
cat( "Generating figure: ", outfile, "\n" )
pdf( file=outfile, width=8, height=5 )
par(mar=c(1,1,1,1))
plot( x=c(0,250000000), y=c(0,-24), ann=F, axes=F, pch=NA )
# for( i in 1:24 ) {
for( i in 1:23 ) {
	# ideogram
	inds <- which(ideo$chr==i)
	for( j in inds ) {
		rect(xleft=ideo$chromStart[j], xright=ideo$chromEnd[j], ybottom=-i+0.4, ytop=-i+0.6, col=gieColor$hex[which(gieColor$name==ideo$gieStain[j])], border=F)
	}
	rect(xleft=min(ideo$chromStart[inds]), xright=max(ideo$chromEnd[inds]), ybottom=-i+0.4, ytop=-i+0.6, border="grey", lwd=0.5)
	# xmap
	inds <- which(xmap$chrom==i)
	for( j in inds ) {
		# rect(xleft=xmap$start[j], xright=xmap$end[j], ybottom=-i+0.0, ytop=-i+0.2, col=rgb(red=105,green=205,blue=205,alpha=90,maxColorValue=255), border=F)	#paleturquoise3
		rect(xleft=xmap$start[j], xright=xmap$end[j], ybottom=-i+0.0, ytop=-i+0.2, col=rgb(red=62,green=158,blue=139,alpha=35/100*255,maxColorValue=255), border=F)	##3E9E8B at 35% opacity
	}
}
mtext( side=2, text=c(1:22,"X"), at=0-(1:23)+0.3, col="paleturquoise4", adj=1, xpd=NA, line=-1, las=1, cex=0.7 )

dev.off()

 

