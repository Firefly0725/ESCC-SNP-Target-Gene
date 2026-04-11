#!/bin/bash

ALL_OVERLAP="SNP_LD.txt"
GENE_BODY="genebody_new1.bed"
PROMOTER_RESULT="type1_LD.txt"
OUTPUT="type2_LD.txt"

get_other_fragments() {
    local bedpe="$1"
    local ref="$2"

    bedtools intersect -a "$bedpe" -b "$ref" -wa -wb |
    cut -f 1-3 | sort -u |
    bedtools intersect -a "$bedpe" -b - -wa -wb |
    cut -f 4-10 | sort -u

    awk '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10}' OFS="\t" "$bedpe" |
    bedtools intersect -a - -b "$ref" -wa -wb |
    cut -f 1-3 | sort -u |
    bedtools intersect -a "$bedpe" -b - -wa -wb |
    cut -f 4-10 | sort -u
}

> "$OUTPUT"

for cell in TT KYSE140 KYSE180 KYSE70; do
    echo "Processing $cell ..."
    BEDPE="${cell}.bedpe"

    bedtools pairtobed -a "$BEDPE" -b "$ALL_OVERLAP" > tmp_${cell}.bedpe

    {
        get_other_fragments tmp_${cell}.bedpe "$ALL_OVERLAP"
    } | sort -u > tmp_frags_${cell}.bed

    bedtools intersect -a tmp_frags_${cell}.bed -b "$GENE_BODY" -F 1.00 -wa -wb |
    cut -f 4-7,12 | sort -u > type2_${cell}_tmp.txt

    grep -Fvxf "$PROMOTER_RESULT" type2_${cell}_tmp.txt > type2_${cell}.txt

    cat type2_${cell}.txt >> "$OUTPUT"
    rm -f tmp_${cell}.bedpe tmp_frags_${cell}.bed type2_${cell}_tmp.txt
done

sort -u "$OUTPUT" -o "$OUTPUT"

echo "Done. Result saved in $OUTPUT"
