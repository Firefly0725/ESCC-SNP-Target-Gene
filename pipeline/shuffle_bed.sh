#!/bin/bash

cell_lines=("KYSE140_H3K27ac" "KYSE70_H3K27ac" "KYSE180_H3K27ac" "TT_H3K27ac" "TE5_H3K27ac" \
"KYSE70_H3K4me1" "TE5_H3K4me1" "KYSE70_H3K4me3" "TE5_H3K4me3" "KYSE140_SOX2" "TE5_SOX2")

for cell in "${cell_lines[@]}"; do
    echo "Processing $cell ..."

    bedpe_file="${cell}.bed"
    snp_file="SNP.txt"
    chrom_sizes="hg19.chrom.sizes"
    out_file="shuffle_${cell}.txt"
    > "$out_file"

    for i in $(seq 1 10000); do
        bedtools shuffle -i "$bedpe_file" -g "$chrom_sizes" -chrom > test.bed
        count=$(intersectBed -a test.bed -b "$snp_file" -wb | cut -f 4-8 | sort -u | wc -l)
        echo "$count" >> "$out_file"
    done
    echo "Finished $cell, results saved to $out_file"
done


cell_lines=("KYSE140_H3K27ac" "KYSE70_H3K27ac" "KYSE180_H3K27ac" "TT_H3K27ac" "TE5_H3K27ac" \
"KYSE70_H3K4me1" "TE5_H3K4me1" "KYSE70_H3K4me3" "TE5_H3K4me3" "KYSE140_SOX2" "TE5_SOX2")

for cell in "${cell_lines[@]}"; do
    echo "Processing $cell ..."

    bedpe_file="${cell}.bed"
    snp_file="SNP_LD.txt"
    chrom_sizes="hg19.chrom.sizes"
    out_file="shuffle_LD_${cell}.txt"
    > "$out_file"

    for i in $(seq 1 10000); do
        bedtools shuffle -i "$bedpe_file" -g "$chrom_sizes" -chrom > test_LD.bed
        count=$(intersectBed -a test_LD.bed -b "$snp_file" -wb | cut -f 4-8 | sort -u | wc -l)
        echo "$count" >> "$out_file"
    done
    echo "Finished $cell, results saved to $out_file"
done
