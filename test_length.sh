#!/bin/bash

#Before
#$1: Query sequence name
#$2: Query sequence length
#$3: Query start coordinate (0-based)
#$4: Query end coordinate (0-based)
#$5: ¡¥+¡¦ if query/target on the same strand; ¡¥-¡¦ if opposite
#$6: Target sequence name
#$8: Target start coordinate on the original strand
#$9: Target end coordinate on the original strand
#$10: Number of matching bases in the mapping
#$11: Number bases, including gaps, in the mapping
#$NF: cigar string

awk '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$NF}' $1 | awk '{split($1,a,"|");$1=a[2];$1=$1"\t"a[3];print $0}' | awk '{gsub(/:/,"\t",$2); print $0;}' | awk '{$6=$6"\t"((($6-$5)/$4)*100); print $0;}' | awk '{$13=$13"\t"(($12/$13)*100); print $0;}' | sort -k 14,14nr -k 7,7nr -k 4,4nr > ../final_paf/$1

#After
#$1: id
#$2: Genus
#$3: Species
#$4: Query sequence length
#$5: Query start coordinate (0-based)
#$6: Query end coordinate (0-based)
#$7: (($5-$4)/$3)*100 query cover
#$8: ¡¥+¡¦ if query/target on the same strand; ¡¥-¡¦ if opposite
#$9: Target sequence name
#$10: Target start coordinate on the original strand
#$11: Target end coordinate on the original strand
#$12: Number of matching bases in the mapping
#$13: Number bases, including gaps, in the mapping
#$14: ($12/$13)*100 identity
#$15: cigar string
