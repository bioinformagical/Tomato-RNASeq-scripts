# README #

This folder contains results from data analysis performed with code in the
parent directory.

What's here:

* CvT-[method]-[assembly].txt - where method is one of "DESeq2" and "edgeR" and assembly is one of "SL4" or "SL5". These were made by Markdown files named FindControVsStressDEGenes-[assembly]-[method].Rmd. These files report differentially expressed genes between treatment and controls, within a single genotype or treatment duration. Genes are from assembly SL4 or SL5.

* CvT-compare-[assembly].txt - where "assembly" is one of "SL4" or "SL5". These files report results from CompareDEGenesByMethod.Rmd, which compares results from using DESeq2 or edgeR to find differentially expressed genes between treatment and controls, within a single genotype or treatment duration. 

* MvW-[assembly].txt and MvW.34.75-[assembly].txt - where "assembly" is one of "SL4" or "SL5." Files named MvW-[assembly].txt report differences in expression between genotypes _are_ and VF36, same temperature and same time point. These files are the output of FindMutantVsWildtypeDEGenes-DESeq2.Rmd.

* muday-144-[assembly]-salmon.txt - where "assembly" is one of "SL4" or "SL5." These files are the output of AddGeneAnnotations-[assembly].Rmd. These report the number RNA-Seq fragments aligned per gene, the "raw" counts, created from input files named muday-144-[assembly]_salmon.merged.gene_counts.tsv, the output of nf-core/rnaseq RNA-Seq data analysis pipeline.

* muday-144-[assembly]-salmon_scaled.txt - where assembly is one of "SL4" or "SL5". These files are the output of InvestigateF3H-and-MakeScaledCounts.Rmd. These files contains scaled gene expression values per gene, expressed in counts per million.

* sample_renaming_summary.txt - output of AddGeneAnnotations-[assembly].Rmd. This file reports how samples were relabeled to correct file name errors introduced by the sequencing provider. 

* MvW-temp-[assembly].xlsx - output of FindTreatmentEffectAcrossGenotypes-DESeq2.Rmd, which finds differentially expressed genes between treatment (28C, 34C) and genotype (are & VF36). The excel files include a different sheet for 
each timepoint (15,30,45,75). The results are ordered by pvalue with the most significant at the top. The method used was DESeq2 and the output includes both assemblies SL4 & SL5.  
