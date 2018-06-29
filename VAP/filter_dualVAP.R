#! /usr/bin/env Rscript

## filter_dualVAP.R
## Eva KF Chan
## https://github.com/ekfchan
## Filters dual-sample VAP output using the following criteria: 
##	INS/DEL:	(Type == "insertion | Type == "deletion) &  Found_in_self_molecules != "no" & Size >= 500 
##	DUP:		Type ~ /dup/ & Found_in_self_molecules != "no" 
## 	INV bkpts:	(Type !~ /_nbase/ & Fail_assembly_chimeric_score != "fail" & Found_in_self_molecules != "no") & Partner == "PASS"
##	TRANS:		Type !~ /_common/ & Type !~ /_segdup/ & Fail_assembly_chimeric_score != "fail" & Found_in_self_molecules <> "no"

## Script is compatible with SMAP File Version 0.7
## Assumes: 
##	1) There is a column header signified by the line beginning with "#h"
##	2) The following column headers are present: Type, Found_in_self_molecules, Size, Fail_assembly_chimeric_score, LinkID, SmapEntryID 

## Arguments: 
## ==========
## Requires: 
## 	input.smap: variants_combine_filters_vs_control_inMoleRefine1.smap from Dual-Sample VAP analysis (Bionano Solve v3.2)
##	output.smap:  Asubset of the input.smap 

args <- commandArgs(TRUE)
if( length( args )!=2 ) { stop("usage: ./filter_dualVAP.R input.smap output.smap\n", call. = FALSE) }

infile = args[1]
outfile = args[2]
cat( "Input:", infile, "\n" )
cat( "Output:", outfile, "\n" )
if(!file.exists(infile)) { stop("Can't file",infile,"\n") }
if(file.exists(outfile)) { stop(outfile,"already exists!\n") }

## Input data: 
## ===========
mycon <- file(infile, open="r")
# Get header first 
while( length(myhead <- readLines(mycon, n=1, warn=F))==1 ) {
	if( substr(myhead,1,2) == "#h" ) {
		break
	}
}
myhead = strsplit(myhead,split="\t")[[1]]
# Slurp in data 
dat <- read.table( file=mycon, header=F, sep="\t", comment.char="#", stringsAsFactors=F, na.strings="", col.names=sub("#h\\s{0,}","",myhead) )
close(mycon)

## Filter Flags:
## =============
## Insertions & Deletions
##	AND( OR(Type = "insertion", Type = "deletion), Found_in_self_molecules <> "no", Size >= 500 )
dat$PASS_insdel = ifelse((dat$Type=="deletion" | dat$Type=="insertion") & dat$Found_in_self_molecules!="no" & dat$Size>=500, TRUE, FALSE)

## Duplications
##	Type ~ /dup/ & Found_in_self_molecules != "no" 
dat$PASS_dup = FALSE
dat$PASS_dup[intersect(grep("dup",dat$type), which(dat$Found_in_self_molecules!="no"))] = TRUE

## Translocations: 
##	Type !~ /_common/ & Type !~ /_segdup/ & Fail_assembly_chimeric_score != "fail" & Found_in_self_molecules <> "no"
dat$PASS_trans = FALSE
dat$PASS_trans[intersect(setdiff(grep("trans",dat$Type),union(grep("common",dat$Type),grep("segdup",dat$Type))), which(dat$Fail_assembly_chimeric_score!="fail" & dat$Found_in_self_molecules!="no"))] = TRUE

## Inversions:
##	breakpoints: (Type !~ /_nbase/ & Fail_assembly_chimeric_score != "fail" & Found_in_self_molecules != "no") & Partner == "PASS"
dat$PASS_invbp = FALSE 
dat$PASS_invbp[intersect(setdiff(grep("inversion",dat$Type),grep("nbase",dat$Type)), which(dat$Fail_assembly_chimeric_score!="fail" & dat$Found_in_self_molecules!="no"))] = TRUE
dat$PASS_inv = dat$PASS_invbp[match(dat$LinkID,dat$SmapEntryID)] & dat$PASS_invbp


## Export filtered: 
## ================
dat$PASS = (dat$PASS_insdel | dat$PASS_dup | dat$PASS_trans | dat$PASS_inv)
write.table(dat[dat$PASS,1:length(myhead)],file=outfile,append=F,quote=F,sep="\t",na="-",row.names=F,col.names=myhead)

q('no')
