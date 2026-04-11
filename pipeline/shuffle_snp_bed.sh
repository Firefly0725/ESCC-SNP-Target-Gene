#!/bin/bash

cell_lines=("KYSE140_H3K27ac" "KYSE70_H3K27ac" "KYSE180_H3K27ac" "TT_H3K27ac" "TE5_H3K27ac" \
"KYSE70_H3K4me1" "TE5_H3K4me1" "KYSE70_H3K4me3" "TE5_H3K4me3" "KYSE140_SOX2" "TE5_SOX2")

for cell in "${cell_lines[@]}"; do
    echo "Processing $cell ..."

    bedpe_file="${cell}.bed"
    snp_file="SNP_overall.txt"
    chrom_sizes="hg19.chrom.sizes"
    out_file="shuffle_snp_${cell}.txt"
    > "$out_file"

    for i in $(seq 1 10000); do
        shuf -n 30 "$snp_file" > random_30.bed
        count=$(intersectBed -a "$bedpe_file" -b random_30.bed -wb | cut -f 4-7 | sort -u | wc -l)
        echo "$count" >> "$out_file"
    done
    echo "Finished $cell, results saved to $out_file"
done

cell_lines=("KYSE140_H3K27ac" "KYSE70_H3K27ac" "KYSE180_H3K27ac" "TT_H3K27ac" "TE5_H3K27ac" \
"KYSE70_H3K4me1" "TE5_H3K4me1" "KYSE70_H3K4me3" "TE5_H3K4me3" "KYSE140_SOX2" "TE5_SOX2")

for cell in "${cell_lines[@]}"; do
    echo "Processing $cell ..."

    bedpe_file="${cell}.bed"
    snp_file="SNP_overall.txt"
    chrom_sizes="hg19.chrom.sizes"
    out_file="shuffle_snp_LD_${cell}.txt"
    > "$out_file"

    for i in $(seq 1 10000); do
        shuf -n 132 "$snp_file" > random_30.bed
        count=$(intersectBed -a "$bedpe_file" -b random_30.bed -wb | cut -f 4-7 | sort -u | wc -l)
        echo "$count" >> "$out_file"
    done
    echo "Finished $cell, results saved to $out_file"
done
