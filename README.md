# SolveMore
Handy scripts for examining Bionano Solve data

_Scripts have been validated for Bionano Solve v3.2 and BNX format v.1.2_

**plotters**

This dir contains scripts for generating plots and figures. 

[**_bionano_ideogram_**](plotters/bionano_ideogram.R)

An R script that takes (a modified) Bionano alignment file (`*.xmap`) as input and generates an _ideogram_ of the data such that the overlaid colour intensity reflects the amount of genome maps aligned to the specific region. It requires a bed file of cytobands to set up the figure background (see a human example in `data/hg38_goldenpath_cytoBandIdeo.txt`). 
The input to the R script is a `bed` format of the `*.xmap`, which can be derived using the command: 
```sh
grep -v "^#" EXP_REFINEFINAL1.xmap | awk 'BEGIN{FS="\t";OFS="\t"}; {print $3, int($6), int($7)}' >  EXP_REFINEFINAL1.xmap_hg38.bed
```
The resulting `bed` can then be used as input to `bionano_ideogram`: 
```R
bionano_ideogram.R hg38_goldenpath_cytoBandIdeo.txt EXP_REFINEFINAL1.xmap_hg38.bed ideogram.pdf
```
Example:  An example is presented as [Fig.1 in the 778 Gen. Res. paper](https://genome.cshlp.org/content/28/5/726.full#F1). 

[**_molecule_diagnostic_plots_**](plottes/molecule_diagnostic_plots.R)

A R script that takes the 0-channel data from a bnx file (v1.3) and generates a series of diagnostic plots on the molecule data. 
The input data can be obtained using a simple `grep`: 
```sh
grep "^0" all_sorted.bnx > all_sorted_0channel.txt
```
The resulting text file can then be used as input to `molecule_diagnostic_plots`: 
```R
Rscript molecule_diagnostic_plots.R all_sorted_0channel.txt all_sorted_0channel_plots.pdf
```
Example:  As example is provided in [`data/molecule_diagnostic_plots.pdf`](data/molecule_diagnostic_plots.pdf)

[**_wgCNV_**](plotters/wgCNV.R)

A R script to plot genome-wide CNV profile, taken as the `fractionalCopyNumber` from `output/contigs/alignmolvref/copynumber/cnv_rcmap_exp.txt` of the Bionano Solve v3 Pipeline. 

_This script is here for historical reasons, as I don't think the plot is very informative...!_



