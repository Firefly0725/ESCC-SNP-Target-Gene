#!/bin/bash

ALL_OVERLAP="SNP_LD.txt"
GENE_BED="gene.bed"
TYPE1_FILE="type1_LD.txt"
TYPE2_FILE="type2_LD.txt"
OUTPUT="type3_LD.txt"

get_other_fragments() {
    local bedpe="$1"
    local ref="$2"

    bedtools intersect -a "$bedpe" -b "$ref" -wa -wb |
    cut -f 1-3 | sort -u |
    bedtools intersect -a "$bedpe" -b - -wa -wb |
    awk '{print $4,$5,$6,$7,$8,$9,$10,$1,$2,$3}' OFS="\t" | sort -u

    awk '{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10}' OFS="\t" "$bedpe" |
    bedtools intersect -a - -b "$ref" -wa -wb |
    cut -f 1-3 | sort -u |
    bedtools intersect -a "$bedpe" -b - -wa -wb |
    awk '{print $4,$5,$6,$7,$8,$9,$10,$1,$2,$3}' OFS="\t" | sort -u
}

> "$OUTPUT"

for cell in TT KYSE140 KYSE180 KYSE70; do
    echo "Processing $cell ..."
    BEDPE="${cell}.bedpe"

    pairToBed -a "$BEDPE" -b "$ALL_OVERLAP" | cut -f 1,3,5,7-10 > test.bed

    bedtools intersect -a test.bed -b "$GENE_BED" -F 1.00 -wa -wb |
    cut -f 4-11 | sort -u > type3.bed

    pairToBed -a "$BEDPE" -b "$ALL_OVERLAP" > test1.bed

    {
        get_other_fragments test1.bed "$ALL_OVERLAP"
    } | sort -u > cons.txt

    awk '{print $4,$5,$6,$7,$1,$2,$3,$8,$9,$10}' OFS="\t" cons.txt > cons.bed

    bedtools intersect -a type3.bed -b cons.bed -wa -wb |
    cut -f 1-8,13-18 | sort -u > type3.txt

    awk '{
        if (($6 > $11 && $7 < $13) || ($6 > $14 && $7 < $10))
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
    }' OFS="\t" type3.txt > type3_1.txt

    awk '{
        if ($3 < $6) print $1,$2,$3,$4,$8,$6-$3
        if ($2 > $7) print $1,$2,$3,$4,$8,$2-$7
    }' OFS="\t" type3_1.txt > type3_2.txt

    cut -f 1-5 type3_2.txt | sort -u | while read key; do
        grep -F "$key" type3_2.txt | sort -n -k6 | head -n1 >> type3_$cell.txt
    done

    rm -f test.bed type3.bed test1.bed cons.txt cons.bed type3.txt type3_1.txt type3_2.txt
done

cat type3_TT.txt type3_KYSE140.txt type3_KYSE70.txt type3_KYSE180.txt |
cut -f 1-5 | sort -u > type3_merged.txt

grep -Fvxf "$TYPE1_FILE" type3_merged.txt > type3_tmp.txt
grep -Fvxf "$TYPE2_FILE" type3_tmp.txt > "$OUTPUT"

rm -f type3_merged.txt type3_tmp.txt

echo "Done. Result saved in $OUTPUT"
