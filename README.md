# SolveMore
Handy scripts for examining Bionano Solve data

_Scripts have been validated for Bionano Solve v3.2, BNX format v.1.2 and SMAP format v0.7_

## **plotters**

This dir contains scripts for generating plots and figures. 

[**_bionano_ideogram_**](plotters/bionano_ideogram.R)

An R script that takes (a modified) Bionano alignment file (`*.xmap`) as input and generates an _ideogram_ of the data such that the overlaid colour intensity reflects the amount of genome maps aligned to the specific region. It requires a bed file of cytobands to set up the figure background (see a human example in `data/hg38_goldenpath_cytoBandIdeo.txt`). 
The input to the R script is a `bed` format of the `*.xmap`, which can be derived using the command: 
```sh
grep -v "^#" EXP_REFINEFINAL1.xmap | awk 'BEGIN{FS="\t";OFS="\t"}; {print $3, int($6), int($7)}' >  EXP_REFINEFINAL1.xmap_hg38.bed
```
The resulting `bed` can then be used as input to `bionano_ideogram`: 
```sh
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
```sh
Rscript molecule_diagnostic_plots.R all_sorted_0channel.txt all_sorted_0channel_plots.pdf
```
Example:  As example is provided in [`data/molecule_diagnostic_plots.pdf`](data/molecule_diagnostic_plots.pdf)

[**_wgCNV_**](plotters/wgCNV.R)

A R script to plot genome-wide CNV profile, taken as the `fractionalCopyNumber` from `output/contigs/alignmolvref/copynumber/cnv_rcmap_exp.txt` of the Bionano Solve v3 Pipeline. 

_This script is here for historical reasons, as I don't think the plot is very informative...!_

## **manipulators**

This dir contains scripts for post-analysis and manipulation of Pipeline outputs. 

[**_tally_SVtypes.sh_**](manipulators/tally_SVtypes.sh)

A shell script to tally up the different SV types from a smap file, printing the two-column results to stdout. It correctly accounts for the double entries for each inversions. It does not discriminate (exclude) `_nbase`, `_common`, or `_segdup` types (filtering scripts can be used for this). 

```sh
tally_SVtypes.sh exp_refineFinal1_merged_filter_inversions.smap
insertions	11068
deletions	9182
duplications	80
inversions	213
intra-chr	23
inter-chr	84
TOTAL SV	20650
```

[**_filter_dualVAP_**](manipulators/filter_dualVAP.R)

A R script to filter the "raw" variants_combine_filters_vs_control_inMoleRefine1.smap file from a Dual-Sample VAP analysis, returing a reduced smap. In brief, SVs of types `_nbase`, `_common`, and `_segdup` are excluded, as are SVs failing Assembly Chimeric Score or not found in self-molecules, and insertions/deletions smaller than 500 bp in size. For inversions, only those where BOTH breakpoints pass inclusion criteria are returned. 

```sh
Rscript filter_dualVAP.R variants_combine_filters_vs_control_inMoleRefine1.smap filteredVAP.smap
```

