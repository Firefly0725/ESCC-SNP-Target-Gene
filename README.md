# ESCC-SNP-Target-Gene
The pipeline includes scripts for defining SNP–target gene interactions and integrating multi-omics datasets as described in the method details.

## 1 Pipeline
### 1.1 Identify SNP-Target Gene Pairs
Based on intra-chromosomal interaction data derived from HiChIP processing and the definitions of the three types of SNP–target gene pairs, we identified potential interaction pairs. This was primarily done using BedTools, along with custom scripts.

Based on the interaction location between SNPs and genes, potential target genes were categorized into three types: Type 1, Type 2, and Type 3. Chromatin interactions between SNPs and genomic fragments were obtained from HiChIP data.  

* **type1.sh**  
  If the promoter region, as defined in the gene structure, is entirely located within the chromatin fragment interacting with the SNP, the regulatory relationship is classified as Type 1. 
* **type2.sh**  
  If the gene body region (excluding the promoter) overlaps with the SNP-interacting fragment, but the promoter does not, the relationship is defined as Type 2.
* **type3.sh**  
  If the entire gene locus lies between the SNP and the interacting fragment, and there is no direct contact between the SNP and the gene, the pair is classified as Type 3.

This categorization was automatically performed using the interaction peaks identified from HiChIP data, based on the above definitions.

### 1.2 DNA Shuffling Experiment
To evaluate whether ESCC-related SNPs are significantly enriched in specific genomic regions or histone-binding regions, we performed two groups of DNA shuffling experiments using the shuffle function from BedTools.

**The first experiment** involved random selection of peak regions from HiChIP and ChIP-seq data. The main steps were as follows: 

(1) Based on the peak regions from ChIP-seq and HiChIP data, as well as the positions of promoters and gene bodies, chromosomal segments of the same length were randomly selected from the same chromosome. 

(2) The number of overlaps between the 30 ESCC-related SNPs and the randomly selected segments from each dataset was calculated. 

(3) This process was repeated 10,000 times, and the number of overlaps from each iteration was recorded. To reduce computational time, only 1,000 iterations were performed for HiChIP data.

**The second experiment** involved random selection of SNPs. The main steps were as follows: 

(1) Based on GWAS results for ESCC from the GWAS Atlas database, 30 SNPs were randomly selected from a pool of 570,000 variant sites. 

(2) The number of overlaps between these randomly selected SNPs and ChIP-seq or HiChIP peak regions was calculated. 

(3) This process was repeated 10,000 times, with overlap counts recorded for each iteration.

For the enrichment significance analysis of LD SNPs, only the first round of experiments was conducted. Unlike the analysis of ESCC SNPs, 132 SNP loci were randomly selected in this case. 

* **shuffle_bedpe.sh / shuffle_bed.sh / shuffle_gene.sh**
* **shuffle_snp_bedpe.sh / shuffle_snp_bed.sh**

### 1.3 Validate SNP-Target Gene Pairs
* **eQTL.R**
  eQTL validation was performed using esophageal tissue data from the GTEx database (v8). The get_qtls function from the R package Qtlizer was employed to batch-retrieve eQTL information for the SNPs used in this study, returning associated genes and significance levels for each SNP across all available tissues. The eQTL results were then inner-joined with the originally identified SNP–target gene pairs based on SNP identifier and gene symbol, retaining only concordant records. The resulting associations were further filtered to include only those containing the keyword "Esophagus Mucosa" and "Esophagus Muscularis" in the tissue field, ensuring that the regulatory evidence was specific to esophageal tissue.

## 2 Figure
