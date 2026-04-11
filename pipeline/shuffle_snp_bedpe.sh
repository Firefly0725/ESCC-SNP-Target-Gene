#!/bin/bash

cell_lines=("KYSE140" "KYSE70" "KYSE180" "TT")

for cell in "${cell_lines[@]}"; do
    echo "Processing $cell ..."

    bedpe_file="${cell}.bedpe"
    snp_file="SNP_overall.txt"
    chrom_sizes="hg19.chrom.sizes"
    out_file="shuffle_snp_${cell}.txt"
    > "$out_file"

    for i in $(seq 1 1000); do
        shuf -n 30 "$snp_file" > random_30_bedpe.bed
        count=$(pairToBed -a "$bedpe_file" -b random_30_bedpe.bed | cut -f 7-10 | sort -u | wc -l)
        echo "$count" >> "$out_file"
    done
    echo "Finished $cell, results saved to $out_file"
done

cell_lines=("KYSE140" "KYSE70" "KYSE180" "TT")

for cell in "${cell_lines[@]}"; do
    echo "Processing $cell ..."

    bedpe_file="${cell}.bedpe"
    snp_file="SNP_overall.txt"
    chrom_sizes="hg19.chrom.sizes"
    out_file="shuffle_snp_LD_${cell}.txt"
    > "$out_file"

    for i in $(seq 1 1000); do
        shuf -n 132 "$snp_file" > random_30_bedpe.bed
        count=$(pairToBed -a "$bedpe_file" -b random_30_bedpe.bed | cut -f 7-10 | sort -u | wc -l)
        echo "$count" >> "$out_file"
    done
    echo "Finished $cell, results saved to $out_file"
done
