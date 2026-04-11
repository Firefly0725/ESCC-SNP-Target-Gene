#!/bin/bash

cell_lines=("KYSE140" "KYSE70" "KYSE180" "TT")

for cell in "${cell_lines[@]}"; do
    echo "Processing $cell ..."

    bedpe_file="${cell}.bedpe"
    snp_file="SNP.txt"
    chrom_sizes="hg19.chrom.sizes"
    out_file="shuffle_${cell}.txt"
    > "$out_file"

    for i in $(seq 1 1000); do
        bedtools shuffle -i "$bedpe_file" -g "$chrom_sizes" -chrom -bedpe > test.bedpe
        count=$(pairToBed -a test.bedpe -b "$snp_file" | cut -f 7-11 | sort -u | wc -l)
        echo "$count" >> "$out_file"
    done
    echo "Finished $cell, results saved to $out_file"
done


cell_lines=("KYSE140" "KYSE70" "KYSE180" "TT")

for cell in "${cell_lines[@]}"; do
    echo "Processing $cell ..."

    bedpe_file="${cell}.bedpe"
    snp_file="SNP_LD.txt"
    chrom_sizes="hg19.chrom.sizes"
    out_file="shuffle_LD_${cell}.txt"
    > "$out_file"

    for i in $(seq 1 1000); do
        bedtools shuffle -i "$bedpe_file" -g "$chrom_sizes" -chrom -bedpe > test_LD.bedpe
        count=$(pairToBed -a test_LD.bedpe -b "$snp_file" | cut -f 7-11 | sort -u | wc -l)
        echo "$count" >> "$out_file"
    done
    echo "Finished $cell, results saved to $out_file"
done
