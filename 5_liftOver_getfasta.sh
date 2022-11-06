#!/bin/bash

## mouse to human

liftOver NCCE_NR2F2_motif_cons_mouse.bed mm10ToHg38.over.chain.gz NCCE_NR2F2_motif_cons_human.bed coup_unmapped

bedtools getfasta -fi hg38.fa -bed NCCE_NR2F2_motif_cons_human.bed -fo NCCE_NR2F2_motif_cons_human.fa


liftOver NCCE_NKX_motif_cons_mouse.bed mm10ToHg38.over.chain.gz NCCE_NKX_motif_cons_human.bed nkx_unmapped




## mouse to zebrafish

liftOver NCCE_NR2F2_motif_cons_mouse.bed mm10ToDanRer11.over.chain.gz NCCE_NR2F2_motif_cons_zebrafish.bed coup_unmapped

bedtools getfasta -fi danRer11.fa -bed NCCE_NR2F2_motif_cons_zebrafish.bed -fo NCCE_NR2F2_motif_cons_zebrafish.fa


liftOver NCCE_NKX_motif_cons_mouse.bed mm10ToDanRer11.over.chain.gz NCCE_NKX_motif_cons_zebrafish.bed nkx_unmapped
