#!/bin/bash

ml biology bedtools

bedtools intersect -a NCCE_NKX_regions_for_cons.bed -b ../../phastConsElements60wayPlacental_score300.bed > NCCE_NKX_mouse_phastConsPlacental_score300.bed
