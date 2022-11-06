library(seqinr)
library(reshape2)

motif = "GGTC(A|G)|(T|C)GACC"

motif.len = 5

fa.in = "phastConsElements60wayPlacental_score300.fa"
bed.out = "phastConsElements60wayPlacental_score300_NR2F2_2motifs_500bp_flank.bed"



fa = read.fasta(fa.in, as.string=T)

fa = unlist(lapply(fa, function(x) x[1]))

hits = lapply(fa, function(x) words.pos(motif, x, ignore.case=T))

names(hits) = paste(paste0("peak", 1:length(hits)), names(hits), sep="-")


hits.table = unique(melt(do.call(rbind, hits))[,-2])   # ignore warning message
colnames(hits.table) = c("peak", "motif")
hits.table = hits.table[order(hits.table$peak),]
    
chr = gsub(".+-(.+):.+", "\\1", hits.table$peak)
start = as.numeric(gsub(".+:(.+)-.+", "\\1", hits.table$peak))
end = as.numeric(gsub(".+-", "", hits.table$peak))

mstart = start + hits.table$motif
mend = start+hits.table$motif + motif.len - 1

mpos = paste(chr, paste(mstart, mend, sep="-"), sep=":")


neighbor.coord = data.frame(chr, start=mstart-250, end=mend+250)

write.table(neighbor.coord, bed.out, sep="\t", col.names=F, row.names=F, quote=F)


neighbor.small = data.frame(chr, start=mstart-50, end=mend+50)

write.table(neighbor.small, "phastConsElements60wayPlacental_score300_NR2F2_2motifs_100bp_flank.bed", sep="\t", col.names=F, row.names=F, quote=F)
