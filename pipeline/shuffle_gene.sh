#!/bin/bash

cell_lines=("promoter_new1" "genebody_new1")

for cell in "${cell_lines[@]}"; do
    echo "Processing $cell ..."

    bedpe_file="${cell}.bed"
    snp_file="SNP.txt"
    chrom_sizes="hg19.chrom.sizes"
    out_file="shuffle_${cell}.txt"
    > "$out_file"

    for i in $(seq 1 10000); do
        bedtools shuffle -i "$bedpe_file" -g "$chrom_sizes" -chrom > test_gene.bed
        count=$(intersectBed -a test_gene.bed -b "$snp_file" -wb | cut -f 6-10 | sort -u | wc -l)
        echo "$count" >> "$out_file"
    done
    echo "Finished $cell, results saved to $out_file"
done


for cell in "${cell_lines[@]}"; do
    echo "Processing $cell ..."

    bedpe_file="${cell}.bed"
    snp_file="SNP_LD.txt"
    chrom_sizes="hg19.chrom.sizes"
    out_file="shuffle_LD_${cell}.txt"
    > "$out_file"

    for i in $(seq 1 10000); do
        bedtools shuffle -i "$bedpe_file" -g "$chrom_sizes" -chrom > test_gene.bed
        count=$(intersectBed -a test_gene.bed -b "$snp_file" -wb | cut -f 6-10 | sort -u | wc -l)
        echo "$count" >> "$out_file"
    done
    echo "Finished $cell, results saved to $out_file"
done
