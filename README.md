# ESCC-SNP-Target-Gene
The pipeline includes scripts for defining SNP–target gene interactions and integrating multi-omics datasets as described in the method details.

## 1 Pipeline
### 1.1 Identify SNP-Target Gene Pairs
Based on intra-chromosomal interaction data derived from HiChIP processing and the definitions of the three types of SNP–target gene pairs, we identified potential interaction pairs. This was primarily done using the pairToBed function from BedTools, along with custom scripts.

Based on the interaction location between SNPs and genes, potential target genes were categorized into three types: Type 1, Type 2, and Type 3. Chromatin interactions between SNPs and genomic fragments were obtained from HiChIP data.  

* **type1.sh**  
  If the promoter region, as defined in the gene structure, is entirely located within the chromatin fragment interacting with the SNP, the regulatory relationship is classified as Type 1. 
* **type2.sh**  
  If the gene body region (excluding the promoter) overlaps with the SNP-interacting fragment, but the promoter does not, the relationship is defined as Type 2.
* **type3.sh**  
  If the entire gene locus lies between the SNP and the interacting fragment, and there is no direct contact between the SNP and the gene, the pair is classified as Type 3.

This categorization was automatically performed using the interaction peaks identified from HiChIP data, based on the above definitions.

### 1.2 DNA Shuffling Experiment
* **shuffle_bedpe.sh / shuffle_snp_bedpe.sh**
* **shuffle_bed.sh / shuffle_snp_bed.sh / shuffle_gene.sh**

### 1.3 Validate SNP-Target Gene Pairs
* **eQTL.R**

## 2 Figure
