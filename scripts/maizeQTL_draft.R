rm(list = ls())
library(GENESPACE)
library(ggplot2)
wd <- "/Users/lovell/Desktop/manuscripts/genespace_2022/GENESPACE_data"
load(file.path(wd, "results", "gparMaize.rda"))
gpar <- gparMaize

wd2 <- file.path(wd, "results", "maize")
gpar$paths$orthogroupsDir <- file.path(wd2, basename(gpar$paths$orthogroupsDir))
gpar$paths$orthofinder <- file.path(wd2, basename(gpar$paths$orthofinder))
gpar$paths$results <- file.path(wd2, basename(gpar$paths$results))
gpar$paths$blastDir <- file.path(wd2, basename(gpar$paths$blastDir))
gpar$params$wd <- wd2

# -- Read in the qtl and reformat
qtl <- fread(file.path(wd, "..", "analysis","Li2016_qtlHits.csv"))
par26 <- colnames(qtl)[which(colnames(qtl) == "B97"):ncol(qtl)]
qtl[,u := paste(Trait, Chromosome, Name, lowerCI)]
qtll <- melt(qtl, measure.vars = par26, variable.name = "genome", value.name = "eff")
setkey(qtll, Trait, Chromosome, Name, lowerCI, eff)

# -- test for upper and lower outliers
qtll[,`:=`(outlp = grubbs.test(eff)$p.value,
           outDir = grubbs.test(eff)$alternative),
     by = "u"]

# -- calculate bonf p-val and subset to significant outlier QTL
qtll[,sig := outlp <= 0.05/nrow(qtl)]
qtlOtl <- subset(qtl, u %in% qtll$u[qtll$sig])
qtlOtl[,genome := c("Ki11", "Mo18W", "Mo18W")]

# get v2 genes in each window
hitsf <- file.path(
  wd, "results", "maize_liftover", "results", "v2_v5_synHits.txt.gz")
gfff <- file.path(
  wd, "results", "maize_liftover", "results", "gffWithOgs.txt.gz")
hits <- fread(hitsf, showProgress = F, na.strings = c("", "NA"))
gff <- fread(gfff, showProgress = F, na.strings = c("", "NA"))

qtlRegl <- lapply(1:nrow(qtlOtl), function(i){
  chri <- qtlOtl$Chromosome[i]
  sti <- qtlOtl$lowerCI[i]*1e6
  eni <- qtlOtl$upperCI[i]*1e6
  return(subset(gff, genome == "v2" & chr == chri & end >= sti & start <= eni))
})
names(qtlRegl) <- qtlOtl$u

# get v5 genes in each window
qtlReglv5 <- lapply(qtlRegl, function(x){
  h <- subset(hits, ofID1 %in% x$ofID & isOg & isAnchor)
  rng <- range(with(gff, which(genome == "v5" & ofID %in% h$ofID2)))
  return(gff[rng[1]:rng[2],])
})

# get qtlIntervals for each
qtlInt <- data.table(
  genome = "B73",
  chr = sapply(qtlReglv5, function(x) x$chr[1]),
  start = sapply(qtlReglv5, function(x) min(x$start)),
  end = sapply(qtlReglv5, function(x) max(x$end)),
  col = c("dodgerblue4","dodgerblue4", "dodgerblue4"))


# pull the b73-anchored pangenome for these regions
pgf <- file.path(
  wd, "results", "maize", "results", "B73_pangenomeDB.txt.gz")
pg <- fread(pgf, showProgress = F, na.strings = c("", "NA"))
gfile <- file.path(
  wd, "results", "maize", "results", "gffWithOgs.txt.gz")
pgff <- fread(gfile, showProgress = F, na.strings = c("", "NA"))
pfile <- file.path(
  wd, "results", "maize", "peptide", "Mo18W.fa")
pepMo18W <- readAAStringSet(filepath = pfile)
setorder(pg, pgChr, pgOrd, na.last = T)
pgRegl <- lapply(qtlReglv5, function(x){
  rng <- range(with(pg, which(id %in% x$id)))
  uids <- pg[rng[1]:rng[2],]
  u <- with(uids, unique(paste(pgChr, pgOrd, pgID)))
  pgo <- subset(pg, paste(pgChr, pgOrd, pgID) %in% u)
  return(pgo)
})

pgGenomes <- gsub("^o", "O", par26)
pgGenomes[pgGenomes == "MS71"] <- "Ms71"
pgGenomes[pgGenomes == "NC305"] <- "NC350"
pgWl <- rbindlist(lapply(1:length(pgRegl), function(i){
  x <- data.table(pgRegl[[i]])
  x <- subset(x, isArrayRep & !isNSOrtho & !isSynOgOnly)
  repGenome <- qtlOtl$genome[i]
  nonRepGenome <- pgGenomes[pgGenomes != repGenome]
  x[,hasOtlGenome := any(genome == repGenome), by = c("pgChr", "pgOrd", "pgID")]
  x[,propGenome := round(uniqueN(genome[genome %in% pgGenomes])/length(pgGenomes),2),
    by = c("pgChr", "pgOrd", "pgID")]
  x[,hasRefGenome := any(genome == "B73"), by = c("pgChr", "pgOrd", "pgID")]
  x <- subset(x, hasRefGenome | hasOtlGenome)
  x[,isPrivateOtl := propGenome < .1 & hasOtlGenome]
  dc <- dcast(
    x,
    pgChr + pgOrd + pgID + propGenome + hasOtlGenome + isPrivateOtl ~ genome,
    value.var = "id",
    fun.aggregate = length)
  m <- melt(dc, measure.vars = unique(x$genome), variable.name = "genome", value.name = "CNV")

  setkey(m, pgOrd)
  m[,pgi := paste(pgID, pgOrd)]
  m[,pgf := as.numeric(factor(pgi, levels = unique(pgi)))]
  tit <- with(qtlOtl, sprintf("%s (%s): %s, %s-%s", Trait, genome, Chromosome, round(lowerCI,1), round(upperCI,1)))
  m[,qtl := tit[i]]
  return(m)
}))
pgoo <- merge(pg, pgWl[,c("qtl", "pgChr", "pgOrd", "pgID")],
              by = c("pgChr", "pgOrd", "pgID"), allow.cartesian = T)
dco <- dcast(pgoo, qtl + pgChr + pgOrd + pgID ~ genome, value.var = "id",
             fun.aggregate = function(x) list(unique(x)))
setkey(dco, pgChr, pgOrd, pgID)
fwrite(
  dco,
  file = "/Users/lovell/Desktop/manuscripts/genespace_2022/siData/siData8_maizeQtlPangenome.txt",
  sep = "\t", quote = F)
qtlOtl$genome <- as.character(qtlOtl$genome)
colnames(qtlOtl) <- gsub("^o", "O",colnames(qtlOtl))
colnames(qtlOtl)[colnames(qtlOtl) == "MS71"] <- "Ms71"
colnames(qtlOtl)[colnames(qtlOtl) == "NC305"] <- "NC350"


pgWl$CNV[pgWl$CNV > 2] <- 2
gids <- gpar$genomes$genomeIDs
pgGenomes <- colnames(qtlOtl)[colnames(qtlOtl) %in% gids]

pgWl[,lev := factor(genome, levels = rev(gids))]
pgWl[,isRefOtl := grepl(genome, qtl), by = c("qtl", "genome")]

pdf(file.path(wd, "raw_plots/Fig2dSourcePlot_qtlPAV.pdf"), height = 2, width = 5)
ggplot(subset(pgWl, qtl != "WWASI (Mo18W): 3, 166.8-170.5"), aes(x = as.factor(pgf), y = lev, fill = as.factor(CNV)))+
  geom_tile(aes(alpha = isRefOtl))+
  scale_fill_manual(values = c("lightgrey", "dodgerblue", "darkblue"), guide = "none")+
  scale_alpha_manual(values = c(.8,1), guide = "none")+
  facet_grid(.~qtl, scale = "free", space = "free")+
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "pangenome order")
dev.off()

pdf(file.path(wd, "raw_plots/FigS5SourcePlot_qtlPAV_lgChr3qtl.pdf"), height = 4, width = 8)
ggplot(subset(pgWl, qtl == "WWASI (Mo18W): 3, 166.8-170.5"), aes(x = as.factor(pgf), y = lev, fill = as.factor(CNV)))+
  geom_tile(aes(alpha = isRefOtl))+
  scale_fill_manual(values = c("lightgrey", "dodgerblue", "darkblue"), guide = "none")+
  scale_alpha_manual(values = c(.8,1), guide = "none")+
  facet_grid(.~qtl, scale = "free", space = "free")+
  theme(panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "pangenome order")
dev.off()


# -- get candidate genes
privCand <- subset(pgWl, isPrivateOtl & isRefOtl)
tmp <- unique(with(privCand, paste(genome, pgID)))
cand <- subset(pg, paste(genome, pgID) %in% tmp)

ki11Func <- fread(
  "/Users/lovell/Desktop/manuscripts/genespace_2022/analysis/Zm-Ki11-REFERENCE-NAM-1.0_Zm00030ab.1.interproscan.tsv.gz",
  sep = "\t", header = F, fill = T)
Mo18WFunc <- fread(
  "/Users/lovell/Desktop/manuscripts/genespace_2022/analysis/Zm-Mo18W-REFERENCE-NAM-1.0_Zm00034ab.1.interproscan.tsv.gz",
  sep = "\t", header = F, fill = T)
ki11Func <- subset(ki11Func, grepl(paste(cand$id, collapse = "|"), V1))
Mo18WFunc <- subset(Mo18WFunc, grepl(paste(cand$id, collapse = "|"), V1))

pfile <- file.path(
  wd, "results", "maize", "peptide", "Mo18W.fa")
pepMo18W <- readAAStringSet(filepath = pfile)
pepMo18W <- pepMo18W[cand$id[cand$id %in% names(pepMo18W)]]
pfile <- file.path(
  wd, "results", "maize", "peptide", "Ki11.fa")
pepKi11 <- readAAStringSet(filepath = pfile)
pepKi11 <- pepKi11[cand$id[cand$id %in% names(pepKi11)]]
writeXStringSet(c(pepKi11, pepMo18W),
                filepath = "/Users/lovell/Desktop/manuscripts/genespace_2022/analysis/combCandPep.fa")

tpm1 <- melt(q1, measure.vars = pgGenomes)[,c("value","variable")]
tpm1 <- rbind(tpm1, data.table(value = 0, variable = "B73"))
tpm1[,lev := factor(variable, levels = levels(tp1$lev))]

bf <- file.path(
  wd, "results", "maize", "results", "syntenicBlocks.txt.gz")

blks <- fread(bf, showProgress = F, na.strings = c("", "NA"))
chr3inv <- subset(blks, gen1 == "B73" & gen2 == "Mo18W" & chr1 == "chr3" &
                    orient == "-" & startBp1 > 170e6)
with(chr3inv, endBp1 - startBp1)/1e6

library(ape)
tre <- read.tree("")

p2 <- ggplot(tpm1, aes(x = value, y = lev))+
  geom_bar(stat = "identity")+
  scale_y_discrete(drop=FALSE)

qtlOtl$genome <- as.character(qtlOtl$genome)
q1 <- qtlOtl[2,]
colnames(q1) <- gsub("^o", "O",colnames(q1))
colnames(q1)[colnames(q1) == "MS71"] <- "Ms71"
colnames(q1)[colnames(q1) == "NC305"] <- "NC350"


pgWl$CNV[pgWl$CNV > 2] <- 2
gids <- unique(as.character(pgWl$genome))
pgGenomes <- colnames(q1)[colnames(q1) %in% gids]
tp1 <- subset(pgWl, grepl("WWASI", qtl))

tp1[,lev := factor(genome, levels = c(
  "B73", unique(gids[!gids %in% c("B73", q1$genome)]), q1$genome))]
p3 <- ggplot(tp1, aes(x = as.factor(pgf), y = lev, fill = CNV))+
  geom_tile()+
  scale_fill_viridis_c(guide = "none")

tpm1 <- melt(q1, measure.vars = pgGenomes)[,c("value","variable")]
tpm1 <- rbind(tpm1, data.table(value = 0, variable = "B73"))
tpm1[,lev := factor(variable, levels = levels(tp1$lev))]

p4 <- ggplot(tpm1, aes(x = value, y = lev))+
  geom_bar(stat = "identity")+
  scale_y_discrete(drop=FALSE)

pdf("raw_plots/Fig2dSourcePlot_pavQTL.pdf", height = 4, width = 8)
grid.arrange(arrangeGrob(p1, p2, ncol=2, nrow=1, widths=c(2,1)),
             arrangeGrob(p3, p4, ncol=2, nrow=1, widths=c(2,1)))
dev.off()

g <- file.path(
  wd, "results", "maize", "results", "gffWithOgs.txt.gz")

g <- fread(g, showProgress = F, na.strings = c("", "NA"))
ki11Genes <- fread("/Users/lovell/Downloads/Zm-Ki11-REFERENCE-NAM-1.0_Zm00030ab.1.interproscan.tsv.gz", sep = "\t", header = F, fill = T)


# # pdf("../../figures/Fig2c_maizerip_rawsourceplot.pdf", height = 6, width = 5)
# plot_riparianHits(
#   gpar,
#   useBlks = F,
#   reorderChrs = F,
#   gapProp = 0.005,
#   blackBg = F,
#   genomeIDs = rev(gpar$genomes$genomeIDs),
#   labelTheseGenomes = "B73",
#   refChrCols = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
#                  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"),
#   chrLabCex = .4, genomeLabCex = .6)
# # dev.off()



pdf(file.path(wd, "raw_plots/Fig2cSourcePlot_ripMaize_2.pdf"), height = 6, width = 5)
plot_riparianHits(
  gpar,
  useBlks = F,
  reorderChrs = F,
  gapProp = 0.005,
  blackBg = F,
  genomeIDs = rev(gpar$genomes$genomeIDs),
  labelTheseGenomes = "B73",
  refChrCols = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                 "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"),
  chrLabCex = .4, genomeLabCex = .6, annotatePlot = F)
plot_riparianHits(
  gpar,
  useBlks = F,
  reorderChrs = F,
  gapProp = 0.005,
  blackBg = F,
  genomeIDs = rev(gpar$genomes$genomeIDs),
  onlyTheseRegions = qtlInt,
  labelTheseGenomes = "B73",
  chrLabCex = .4, genomeLabCex = .6, add2plot = T)
dev.off()


pdf(file.path(wd, "raw_plots/Fig2dSourcePlot_ripZoomChr3.pdf"), height = 2, width = 5)
plot_riparianHits(
  gpar,
  useBlks = T,
  useOrder = F,
  reorderChrs = F,
  gapProp = 0.005,
  blackBg = F,
  genomeIDs = c("Tzi8","Mo18W", "B73"),
  onlyTheseRegions = data.table(
    genome = "B73", chr = "chr3", start = 0, end = 1e10, col = "#B2DF8A"),
  chrLabCex = .4, genomeLabCex = .6,
  excludeNoRegChr = T, annotatePlot = F)
plot_riparianHits(
  gpar,
  useBlks = T,
  useOrder = F,
  reorderChrs = F,
  gapProp = 0.005,
  blackBg = F,
  genomeIDs = c("Tzi8","Mo18W", "B73"),
  onlyTheseRegions = qtlInt[2:3,],
  chrLabCex = .4, genomeLabCex = .6,
  excludeNoRegChr = T, add2plot = T)
dev.off()

pdf(file.path(wd, "raw_plots/Fig2dSourcePlot_ripZoomChr6.pdf"), height = 2, width = 5)
plot_riparianHits(
  gpar,
  useBlks = T,
  useOrder = F,
  reorderChrs = F,
  gapProp = 0.005,
  blackBg = F,
  genomeIDs = c("Tzi8","Ki11", "B73"),
  onlyTheseRegions = data.table(
    genome = "B73", chr = "chr6", start = 0, end = 1e10, col = "#B2DF8A"),
  chrLabCex = .4, genomeLabCex = .6,
  excludeNoRegChr = T, annotatePlot = F)
plot_riparianHits(
  gpar,
  useBlks = T,
  useOrder = F,
  reorderChrs = F,
  gapProp = 0.005,
  blackBg = F,
  genomeIDs = c("Tzi8","Ki11", "B73"),
  onlyTheseRegions = qtlInt[1,],
  chrLabCex = .4, genomeLabCex = .6,
  excludeNoRegChr = T, add2plot = T)
dev.off()

pdf(file.path(wd, "raw_plots/Fig2dSourcePlot_ripZoomChr6.pdf"), height = 2, width = 5)
plot_riparianHits(
  gpar,
  useBlks = T,
  useOrder = F,
  reorderChrs = F,
  gapProp = 0.005,
  blackBg = F,
  genomeIDs = c("Tzi8","Ki11", "B73"),
  onlyTheseRegions = data.table(
    genome = "B73", chr = "chr6", start = 0, end = 1e10, col = "#B2DF8A"),
  chrLabCex = .4, genomeLabCex = .6,
  excludeNoRegChr = T, annotatePlot = F)
plot_riparianHits(
  gpar,
  useBlks = T,
  useOrder = F,
  reorderChrs = F,
  gapProp = 0.005,
  blackBg = F,
  genomeIDs = c("Tzi8","Ki11", "B73"),
  onlyTheseRegions = qtlInt[1,],
  chrLabCex = .4, genomeLabCex = .6,
  excludeNoRegChr = T, add2plot = T)
dev.off()


pdf(file.path(wd, "raw_plots/FigS5SourcePlot_fullMaizeRip.pdf"), height = 11, width = 10)
plot_riparianHits(
  gpar,
  useBlks = T,
  useOrder = F,
  reorderChrs = F,
  gapProp = 0.01,
  blackBg = F,
  labelTheseGenomes = "B73",
  genomeIDs = rev(gpar$genomes$genomeIDs),
  chrLabCex = .4, genomeLabCex = .6)
dev.off()
