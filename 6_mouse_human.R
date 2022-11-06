## lift mouse bed to human and zebrafish using liftOver (mm10ToHg38, mm10ToDanRer11)

library(seqinr)
library(reshape2)


motif = "ggtca|ggtcg|tgacc|cgacc"
#motif = "gtacg|gtacc|cgtac|ggtac"
#motif = "agtcg|agtcc|cgact|ggact"
#motif = "tggac|tggat|gtcca|atcca"


nr2f2.human = read.table("NCCE_NR2F2_motif_cons_human.bed", sep="\t", header=F)

nr2f2.human$coords = paste(nr2f2.human$V1, paste(nr2f2.human$V2, nr2f2.human$V3, sep="-"), sep=":")
nr2f2.human$diff = nr2f2.human$V3 - nr2f2.human$V2


## keep only NR2F2 sites with correct sequences
fa = read.fasta("NCCE_NR2F2_motif_cons_human.fa", as.string=T)

fa = unlist(lapply(fa, function(x) x[1]))
fa = data.frame(coords = names(fa), fa)

identical(nr2f2.human$coords, fa$coords)


nr2f2.seq = data.frame(nr2f2.human, fa)


nr2f2.match = nr2f2.seq[grep(motif, nr2f2.seq$fa),]
nr2f2.match = nr2f2.match[which(nr2f2.match$diff <= 10),]   # drop in frequency after 10 (12 in original analysis)


motif.pos = lapply(nr2f2.match$fa, function(x) words.pos(motif, x, ignore.case=T))
names(motif.pos) = nr2f2.match$coords

motif.pos = unique(melt(do.call(rbind, motif.pos))[,-2])
colnames(motif.pos) = c("coords", "motif.pos")

nr2f2.match = merge(nr2f2.match, motif.pos, by="coords")


## generate exact coords of nr2f2 motif in human
nr2f2.exact = data.frame(chr=nr2f2.match$V1, start=nr2f2.match$V2 + nr2f2.match$motif.pos, end=nr2f2.match$V2 + nr2f2.match$motif.pos + 4)
nr2f2.exact = data.frame(NR2F2.exact = paste(nr2f2.exact$chr, paste(nr2f2.exact$start, nr2f2.exact$end, sep="-"), sep=":"), NR2F2.Region.1 = nr2f2.match$coords)

write.csv(nr2f2.exact, "NCCE_NR2F2_exact_coords_human.csv")


## extend to 550bp
nr2f2.550bp = data.frame(chr=nr2f2.match$V1, start=nr2f2.match$V2-275, end=nr2f2.match$V3+275)

write.table(nr2f2.550bp, "NR2F2_cons_human_550bp_flank.bed", sep="\t", col.names=F, row.names=F, quote=F)


## run homer and generate human table, format table using 3_format_homer_result.r
## human only: filter out NKX sites that are not "syntenic" in mouse

homer.hu = read.csv("NCCE_human_all.csv", row.names=1)
cons.mo = read.csv("NCCE_mouse_conserved.csv")
nr2f2.hu = read.csv("NCCE_NR2F2_exact_coords_human.csv", row.names=1)

homer.hu = merge(homer.hu, nr2f2.hu, by="NR2F2.Region.1")


nr2f2.hu = read.table("NCCE_NR2F2_motif_cons_human.bed", sep="\t", header=F)
nr2f2.hu$nr2f2.hu.coords = paste(nr2f2.hu$V1, paste(nr2f2.hu$V2, nr2f2.hu$V3, sep="-"), sep=":")
nr2f2.hu = nr2f2.hu[,4:5]

hu.nr2f2.mo = merge(homer.hu, nr2f2.hu, by.x="NR2F2.Region.1", by.y="nr2f2.hu.coords")


nkx.hu = read.table("NCCE_NKX_motif_cons_human.bed", sep="\t", header=F)
nkx.hu$nkx.hu.coords = paste(nkx.hu$V1, paste(nkx.hu$V2, nkx.hu$V3, sep="-"), sep=":")

nkx.hu.0 = data.frame(V4 = nkx.hu$V4, nkx.hu.coords = paste(nkx.hu$V1, paste(nkx.hu$V2, nkx.hu$V3-2, sep="-"), sep=":"))
nkx.hu.1 = data.frame(V4 = nkx.hu$V4, nkx.hu.coords = paste(nkx.hu$V1, paste(nkx.hu$V2+1, nkx.hu$V3-1, sep="-"), sep=":"))
nkx.hu.2 = data.frame(V4 = nkx.hu$V4, nkx.hu.coords = paste(nkx.hu$V1, paste(nkx.hu$V2+2, nkx.hu$V3, sep="-"), sep=":"))
nkx.hu.all = rbind(nkx.hu.0, nkx.hu.1, nkx.hu.2)

hu.ncce.mo = merge(hu.nr2f2.mo, nkx.hu.all, by.x="NKX.Region", by.y="nkx.hu.coords")

hu.ncce.mo$NKX.NR2F2.mouse = paste(hu.ncce.mo$V4.y, hu.ncce.mo$V4.x, sep="_")

hu.ncce.mo$NCCE.Gap.true = abs((as.numeric(gsub(".+:(\\d+)-.+", "\\1", hu.ncce.mo$NR2F2.exact)) + 2) - (as.numeric(gsub(".+:(\\d+)-.+", "\\1", hu.ncce.mo$NKX.Region)) + 3))



cons.mo$NR2F2.Region = paste0("mouse_NR2F2_", cons.mo$NR2F2.Region)
cons.mo$NKX.Region = paste0("mouse_NKX_", cons.mo$NKX.Region)

cons.mo$NKX.NR2F2.mouse = paste(cons.mo$NKX.Region, cons.mo$NR2F2.Region, sep="_")



ncce.syn = merge(hu.ncce.mo, cons.mo, by="NKX.NR2F2.mouse")

ncce.syn = ncce.syn[,-c(1,3:6,9,10,15,16,18:21,29)]

colnames(ncce.syn) = c("NKX.Region.hu", "Distance.to.TSS.hu", "Gene.Name.hu", "NCCE.Region.hu", "NKX.Motif.hu", "NKX.Strand.hu", "NR2F2.Region.hu", "NCCE.Gap.hu",  "Distance.to.TSS.ms", "Gene.Name.ms", "NR2F2.Region.ms", "NKX.Region.ms", "NCCE.Gap.ms", "NCCE.Region.ms", "NKX.Motif.ms", "NKX.Strand.ms")

ncce.syn = unique(ncce.syn)

ncce.syn = ncce.syn[,c(14,4,11,7,12,1,9,2,10,3,16,6,13,8,15,5)]

ncce.syn$NCCE.Gap.Equal = apply(ncce.syn, 1, function(x) as.numeric(x[13]) == as.numeric(x[14]))
ncce.syn$NKX.Hamming.Dist = apply(ncce.syn, 1, function(x) sum(unlist(strsplit(x[15], split="")) != unlist(strsplit(x[16], split=""))))

ncce.syn$NKX.Motif.Type[which(nchar(ncce.syn$NKX.Motif.ms)==7)] = "AAGTG"
ncce.syn$NKX.Motif.Type[which(nchar(ncce.syn$NKX.Motif.ms)==8)] = "ATTA"

write.csv(ncce.syn, "NCCE_human_mouse_syntenic_gap250.csv")



