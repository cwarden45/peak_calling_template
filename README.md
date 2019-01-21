### Acknowledgements ###

Code developed when providing analysis support from the City of Hope Integrative Genomics Core with requests / suggestions from *Yapeng Su* (Caltech; Broad/Narrow ChIP-Seq, ATAC-Seq).

I would also like to thank the following users that gave us permission to acknowledge these templates (or scripts similar to these templates) were used for analysis of their data: *Zuoming Sun* (ChIP-Seq), *John Chan* (ATAC-Seq QC), and *Jeremy Jones* (Ben Copeland’s project; narrow ChIP-Seq re-analysis / QC).  There is **1** other internal collaborator that currently did not to be explicitly acknowledged here (at this time). 

In general, some projects have been analyzed or re-analyzed with different staff members (sometimes in different labs/institutes, sometimes within the IGC).  So, the use of data with modified versions of these templates may not be what is eventually used in the associated publication.


### \~Temporary Note\~ ###
**I apologize for the confusion, but I would like to emphasize that these are called “templates” because I almost always have to modify the code for each project (beyond the parameter files), meaning they will be more difficult for other people to use in the same way.**  This was unfortunately not immediately clear to me when I created the templates.

I also believe that the process of writing the scripts for analysis (such as the templates) is very important for the learning process, and it is very important that you understand all the steps for analysis before presenting them in a paper.

So, I will post an update when more specific guidance / suggestions can be provided.  **In the meantime, I apologize, but I cannot provide user support for the templates.**

### Reference Annotations ###

To quantify counts around TSS without peak calling run `create_promoter_gtf.R`.  This creates annotation in a format **[gene symbol]_chr:start-stop:strand**.

### Order to Run Scripts ###

1) `cluster_Bowtie2_alignment.py`

2) `coverage_statistics.py`

3) `create_TDFs.py`

4) `cluster_MACS2.py`

5) `peak_statistics.R`

After merged peak GTF is created from **peak_statistics.R**, you can optionally run **run_ChIPseeker_QC.R**.

6) `cluster_HTseq_counts.py`

Must be run after **peak_statistics.R** to have merged peak GTF.

7) `reformat_HTseq_peak.R` or `reformat_HTseq_TSS.py`

NOTE: If normalizing to canonical non-chrM chromosomes (with 'aligned' CPM), then you also need to run **aligned_read_counts.R** prior to creating count table

**reformat_HTseq_peak.R** adds ChIPseeker nearest gene annotations

8) `QC_FPKM_peak.R` or `QC_FPKM_TSS.R`

9) `differential_binding_analysis_peak.R` or `differential_binding_analysis_TSS.R`


### Dependencies (some optional) ###

*Alignment*

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

*Read Counting*

HTseq: http://www-huber.embl.de/HTSeq/doc/install.html#install

*Feature Annotation*

ChIPseeker: http://bioconductor.org/packages/release/bioc/html/ChIPseeker.html

GenomicAlignments: https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html

*Differential Expression*

edgeR: https://bioconductor.org/packages/release/bioc/html/edgeR.html

limma-voom: https://bioconductor.org/packages/release/bioc/html/limma.html

DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html

qvalue: https://bioconductor.org/packages/release/bioc/html/qvalue.html

*Visualization*

gplots: https://cran.r-project.org/web/packages/gplots/index.html

RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html

heatmap.3: https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R

heatmap.3 example: https://www.biostars.org/p/18211/

IGV / igvtools: http://software.broadinstitute.org/software/igv/

*Functional Enrichment / Motif-Finding*

i-cisTarget: https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/index.php

HOMER: http://homer.ucsd.edu/homer/introduction/programs.html

MEME Suite: http://meme-suite.org/index.html

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name|Name of differential binding comparison (used to name output file)|
|plot_groups|Names of columns in *sample_description_file* to be plotted in QC and differential binding plots.  Use commas to plot multiple groups|
|deg_groups|Names of columns in *sample_description_file* to be plotted in QC and differential binding plots.  Use commas to include multiple variables (for multivariate model or gene list filtering)|
|treatment_group|Treatment group for primary variable; enter *continuous* for a continuous variable and a correlation will be provided instead of a fold-change value.|
|Result_Folder|Path to output folder for selected, final results|
|Alignment_Folder|Path to Bowtie2 Alignments|
|Reads_Folder|Path to Reads for Bowtie2 Alignment|
|Cluster_Email|If running alignment on a cluster, e-mail for notifications|
|dba_pvalue_method|Method to Calculate P-value for Counts in Differenital Binding Analysis.  Can be *edgeR*, *limma-voom*, *DESeq2*, *lm* (linear regression), or *aov* (ANOVA)|
|dba_fdr_method|Method to Calculate FDR for Counts in Differenital Binding Analysis.  Can be *BH* (Benjamini and Hochberg),*q-value*, or *q-lfdr*|
|genome|Name of genome build|
|MEM_Limit|Memory allocated to java or job on cluster|
|Bowtie_Ref| Path to Bowtie ref|
|peakType|Are you calling 'narrow' or 'broad' peaks? Can be 'yes' or 'no'.|
|Threads|Number of Threads for Bowtie2 Alignment|
|PE_Reads|Are you using paired-end reads? Can be 'yes' or 'no'.|
|Remove_Duplicates|Use Picard to remove duplicates? Can be 'yes' or 'no'.|
|promoterGTF|Path to promoter .gtf file for htseq-count; produced using `create_promoter_gtf.R`|
|flankTSS|When creating promoter .gtf, flanking distance from TSS (total sizs will be twice this, for upstream and downstream flanking sequence)|
|mergedPeakGTF|Merged peak from separate peak-calling results, using for htseq-count|
|gtfIF|Tag in .gtf file (peak-based or TSS-based) used for quantification by htseq-count|
|sample_description_file|Name of Sample Description File|
|aligned_stats_file|If calculating aligned FPKM / CPM, file to save canonical chromosome aligned read counts (omitting those specified with **omitChr**)|
|chr_length_file|Table of chromosome lengths|
|omitChr|Chromosomes to omit from aligned read count (such as chrM)|
|FPKM_norm|	How to count number of aligned reads: **aligned**, **quantified**, or **TMM**|
|fpkm_binding_cutoff|Rounding value for FPKM|
|minimum_fraction_bound|Minimum fraction of samples with binding above **fpkm_binding_cutoff**. Filter for differential binding analysis.|
|fold_change_cutoff|Minimum fold-change difference for differential binding of peak/promoter|
|cor_cutoff|If using a continuous variable, minimum absolute correlation for differenital binding of peak/promoter|
|sec_cor_cutoff|If comparing two gene lists and using a continuous variable, minimum absolute correlation for differenital binding of peak/promoter|
|pvalue_cutoff|Maximum p-value for differenital binding of peak/promoter|
|fdr_cutoff|Maximum FDR for differenital binding of peak/promoter|
|sec_fold_change_cutoff|If comparing two gene lists, fold-change threshold for list you want to filter out|
|sec_pvalue_cutoff|If comparing two gene lists, p-value threshold for list you want to filter out|
|sec_fdr_cutoff|If comparing two gene lists, FDR threshold for list you want to filter out|
|secondary_trt|If comparing two gene lists, this is treatment group for the list that you want to filter out; enter *continuous* for a continuous variable and a correlation will be provided instead of a fold-change value (also converts second variable from factor to numeric, even if interaction is set to *no*)|
|interaction| Method for comparing an interaction of two variables.  Can be *model* or *no*|
|fpkm_file_peak|Name of File to Contain log2(FPKM + *rpkm_expression_cutoff*) Values for Merged Peaks|
|CPM_file_peak|Name of File to Contain CPM Values for Merged Peaks|
|counts_file_peak|Name of File to Contain Counts for Merged Peaks|
|fpkm_file_TSS|Name of File to Contain log2(FPKM + *rpkm_expression_cutoff*) Values for TSS+/-*flankTSS*|
|CPM_file_TSS|Name of File to Contain CPM Values for TSS+/-*flankTSS*|
|counts_file_TSS|Name of File to Contain Counts for TSS+/-*flankTSS*|
|cluster_distance| Distance metric for dendrogram.  Can be *Euclidean* or *Pearson_Dissimilarity*|
