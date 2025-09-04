#!/bin/bash

## Script detecting miRNA binding sites in input sequences

## input:
# circRNA sequences file: "/data/backsplice_sequence_1.fa"
# miRNA database: "/data/input/.mature_"$SPECIES".fa"
# miRanda parameters: $SCORE, $ENERGY
#	SCORE: best predictions are obtained with higher score
#	ENERGY: best predictions are obtained with lower energy
# AGO2 binding file (optional): "/data/input/AGO2_binding_sites.bed"
# file of region positions for each circRNA (optional): "/data/sequence_extraction/region_to_extract_1.bed"

## command: /scripts/script_miRNA_detection.sh $MIRNA_FILE $SCORE $ENERGY



MIRNA=$1
SEQUENCE="region_to_extract_1.bed"

cat backsplice_sequence_1.fa | grep ">" > headers.txt

cat backsplice_sequence_1.fa | grep -v ">" | cut -c1-20 > first_20_characters.txt
cat backsplice_sequence_1.fa | grep -v ">" > sequences.txt

paste sequences.txt first_20_characters.txt | sed -e 's/\t//' > sequences_plus_first_20_characters.txt

paste headers.txt sequences_plus_first_20_characters.txt | sed -e 's/\t/\n/' > backsplice_sequence_per_miRNA.fa


# miRanda

miranda ../hsa_mature.fa backsplice_sequence_per_miRNA.fa -sc 140 -en -20 -quiet -out output_miRanda.txt

echo "#circ_id miRNA_id score energy circ_start circ_end miRNA_start miRNA_end align_length circ_align_perc miRNA_align_perc" | sed -e 's/ /\t/g' > output_miRanda_b.txt
cat output_miRanda.txt | grep ">" | grep -v ">>" | sed -e 's/>//' |  awk '{print $2,$1,$3,$4,$7,$8,$5,$6,$9,$11,$10}' | sed -e 's/ /\t/g' >> output_miRanda_b.txt
cat output_miRanda_b.txt | sed -e 's/ +/\t/g; s/#//' > output_miRanda_b_per_R.txt


cat output_miRanda_b_per_R.txt | head -n1 > output_miRanda_c_per_R.txt

#head -n1 output_miRanda_b_per_R.txt > output_miRanda_c_per_R.txt
cat backsplice_circRNA_length_1.txt | grep -v "circ_id" > circRNA_length.txt

for CIRC in $( cat circRNA_length.txt | cut -f1 )
do
	L=$( cat circRNA_length.txt | grep $CIRC | cut -f2 )
	cat output_miRanda_b_per_R.txt | grep $CIRC | awk '$5 <= '$L' {print}' >> output_miRanda_c_per_R.txt
done

mv output_miRanda_c_per_R.txt output_miRanda_per_R.txt

