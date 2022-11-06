#!/bin/bash

in1="phastConsElements60wayPlacental_score300_NR2F2_2motifs_500bp_flank.bed"
out1="phastConsPlacental_score300_mouse.txt"

in2="NR2F2_cons_human_550bp_flank.bed"
out2="NCCE_NKX_human_550bp.txt"

in3="NR2F2_cons_zebrafish_550bp_flank.bed"
out3="NCCE_NKX_zebrafish_550bp.txt"




# mouse
annotatePeaks.pl $in1 mm10 -m NKX_STAMP_AAGTG_base4.37.motif NKX_STAMP_ATTA_base4.37.motif > $out1


#human
annotatePeaks.pl $in2 hg38 -m NKX_STAMP_AAGTG_base4.37.motif NKX_STAMP_ATTA_base4.37.motif > $out2


#zebrafish
annotatePeaks.pl $in3 danRer11 -m NKX_STAMP_AAGTG_base4.37.motif NKX_STAMP_ATTA_base4.37.motif > $out3

