library(stringr)
library(data.table)


prefix = "phastConsPlacental_score300_mouse"
#prefix = "NCCE_NKX_human_550bp"
#prefix = "NCCE_NKX_zebrafish_550bp"

species = "mouse"
#species = "human"
#species = "zebrafish"


motifs = fread(paste0(prefix,".txt"))

motifs = as.data.frame(motifs)

colnames(motifs)[22:23] = c("AAGTG.motif", "ATTA.motif")

out = motifs[-which(motifs$AAGTG.motif=="" & motifs$ATTA.motif==""),]
out$motifs.all = paste(out$AAGTG.motif, out$ATTA.motif, sep=",")


out = out[,c(2:4,10,16,24)]

out = unique(out)
nmotif = str_count(out$motifs.all, "\\(")
table(nmotif)


out$rowid = paste0(paste0("R", seq(1, nrow(out))), "X")

out1 = strsplit(out$motifs.all, split="),")
names(out1) = out$rowid

out1 = unlist(out1)
out1 = data.frame(rowid = names(out1), out1)
out1$rowid = gsub("\\d+$", "", out1$rowid)

motif.all = merge(out, out1, by="rowid")

motif.all$dist = as.numeric(gsub("\\(.+", "", gsub("^,", "", motif.all$out1)))
motif.all$motif = gsub(".+\\(([A-Z]+).+", "\\1", motif.all$out1)


motif.all$NKX.Start = motif.all$Start + motif.all$dist
motif.all$NKX.End = motif.all$Start + motif.all$dist + nchar(motif.all$motif) - 1


if (species == "mouse") {
    motif.all$NR2F2.Start = motif.all$Start + 249
    motif.all$NR2F2.End = motif.all$Start + 253
    motif.all$NCCE.Gap = abs((motif.all$Start + 251) - (motif.all$Start + motif.all$dist + 3))
} else if (species %in% c("human", "zebrafish")) {
    motif.all$NR2F2.Start = motif.all$Start + 275
    motif.all$NR2F2.End = motif.all$End - 275
    motif.all$NCCE.Gap = abs((motif.all$Start + motif.all$End - 1)/2 - (motif.all$Start + motif.all$dist + 3))
}


motif.all$NCCE.Start = apply(motif.all, 1, function(x) min(x[11:14]))
motif.all$NCCE.End = apply(motif.all, 1, function(x) max(x[11:14]))

motif.all$NCCE.Region = paste(motif.all$Chr, paste(motif.all$NCCE.Start, motif.all$NCCE.End, sep="-"), sep=":")
motif.all$NKX.Region = paste(motif.all$Chr, paste(motif.all$NKX.Start, motif.all$NKX.End, sep="-"), sep=":")
motif.all$NR2F2.Region = paste(motif.all$Chr, paste(motif.all$NR2F2.Start, motif.all$NR2F2.End, sep="-"), sep=":")

motif.all$NR2F2.Region.1 = paste(motif.all$Chr, paste(motif.all$NR2F2.Start-1, motif.all$NR2F2.End, sep="-"), sep=":")
## for synteny validation

motif.all$NKX.Strand = gsub(".+,(.),.+", "\\1", motif.all$out1)


## human/zebrafish only

motif.out = motif.all[,c(2,16,17,5,6,20,19,15,18,10,21,22)]

colnames(motif.out)[10] = "NKX.Motif"
motif.out$NCCE.Region = gsub(" ", "", motif.out$NCCE.Region)

write.csv(motif.out, "NCCE_zebrafish_all.csv")
## human/zebrafish: before filtering out non-syntenic NKX motifs (but NR2F2 motifs are of the correct sequences)


## mouse only: intersect NKX motifs with phastCons file using bedtools

nkx.bed = data.frame(chr=motif.all$Chr, start=motif.all$NKX.Start, end=motif.all$NKX.End)

write.table(nkx.bed, "NCCE_NKX_regions_for_cons.bed", col.names=F, row.names=F, quote=F, sep="\t")


nkx.cons = read.table("NCCE_NKX_mouse_phastConsPlacental_score300.bed", header=F, sep="\t")
nkx.cons$coords = paste(nkx.cons$V1, paste(nkx.cons$V2, nkx.cons$V3, sep="-"), sep=":")

cons = motif.all[which(motif.all$NKX.Region %in% nkx.cons$coords),]


motif.out = cons[,c(2,16,17,5,6,20,19,15,18,10,21,22)]

colnames(motif.out)[10] = "NKX.Motif"
motif.out$NCCE.Region = gsub(" ", "", motif.out$NCCE.Region)

write.csv(motif.out, "NCCE_mouse_conserved.csv")



nkx = data.frame(chr=cons$Chr, start=cons$NKX.Start-1, end=cons$NKX.End+1, mouse.nkx=paste0("mouse_NKX_", cons$NKX.Region))
nr2f2 = data.frame(chr=cons$Chr, start=cons$NR2F2.Start-1, end=cons$NR2F2.End+1, mouse.nr2f2=paste0("mouse_NR2F2_", cons$NR2F2.Region))


write.table(nkx, "NCCE_NKX_motif_cons_mouse.bed", col.names=F, row.names=F, sep="\t", quote=F)
write.table(nr2f2, "NCCE_NR2F2_motif_cons_mouse.bed", col.names=F, row.names=F, sep="\t", quote=F)


