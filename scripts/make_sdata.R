################################################################################
# vertebrates, Sdata 1-3

wd <- "/Users/lovell/Desktop/manuscripts/genespace_2022/GENESPACE_data/results"
library(GENESPACE)
library(ggplot2)
load(file.path(wd, "gparVerts.rda"))
gpar <- gparVerts

bn <- file.path(wd, "vertebrates")
gpar$paths$results <- file.path(bn, "results")

gpar$paths$mcscanxCall <- "/Users/lovell/Desktop/programs/MCScanX/"

gpar$paths$blastDir <- file.path(bn, "/orthofinder/Results_Feb12/WorkingDirectory")

gpar$paths$orthogroupsDir <- file.path(bn, "/orthofinder/Results_Feb12/Orthogroups")


pltOrd <- gpar$genomes$genomeIDs[c(1:2, 5:4, 7:6, 3, 8, 11:12, 10:9, 15:17, 13:14)]
gff <- fread(file.path(gpar$paths$results, "gffWithOgs.txt.gz"), na.strings = c("", "NA"), showProgress = F)

blks <- fread(file.path(gpar$paths$results, "syntenicBlocks.txt.gz"), na.strings = c("", "NA"), showProgress = F)


## Make SI data files
convert_pgdb2wide <- function(path2PgDB){
  pg <- fread(path2PgDB, na.strings = c("", "NA"), showProgress = F)
  pg <- subset(pg, !is.na(pgID))
  # -- flag non-syntenic orthogroups
  wh <- which(pg$isNSOrtho)
  pg$id[wh] <- paste0(pg$id[wh], "*")

  # -- flag array reps
  wh <- which(!pg$isArrayRep & !pg$isNSOrtho)
  pg$id[wh] <- paste0(pg$id[wh], "+")

  # -- reshape to wide format
  pgw <- dcast(
    pg,
    pgID + pgChr + pgOrd ~ genome,
    value.var = "id",
    fun.aggregate = function(x) list(x))

  # -- order by pg position
  setorder(pgw, pgChr, pgOrd, na.last = T)
  return(pgw)
}

humanPg <- convert_pgdb2wide(
  file.path(gpar$paths$results, "human_pangenomeDB.txt.gz"))
fwrite(
  humanPg,
  file = "/Users/lovell/Desktop/manuscripts/genespace_2022/siData/siData1_humanAnchoredPangenomeAnnotation.txt",
  sep = "\t", quote = F)

chickPg <- convert_pgdb2wide(
  file.path(gpar$paths$results, "chicken_pangenomeDB.txt.gz"))
fwrite(
  chickPg,
  file = "/Users/lovell/Desktop/manuscripts/genespace_2022/siData/siData2_chickenAnchoredPangenomeAnnotation.txt",
  sep = "\t", quote = F)

bout <- with(blks, data.table(
  genome1 = gen1, genome2 = gen2, chr1 = chr1, chr2 = chr2,
  start1 = startBp1, end1 = endBp1, start2 = startBp2, end2 = endBp2,
  orient = orient, nHits = nHits1))
fwrite(
  bout,
  file = "/Users/lovell/Desktop/manuscripts/genespace_2022/siData/siData3_vertSynBlkCoords.txt",
  sep = "\t", quote = F)

################################################################################
# -- grasses, Si data 4-5
blks <- fread("/Users/lovell/Desktop/manuscripts/genespace_2022/GENESPACE_data/results/grasses/results/syntenicBlocks.txt.gz", na.strings = c("", "NA"), showProgress = F)
bout <- with(blks, data.table(
  genome1 = gen1, genome2 = gen2, chr1 = chr1, chr2 = chr2,
  start1 = startBp1, end1 = endBp1, start2 = startBp2, end2 = endBp2,
  orient = orient, nHits = nHits1))
fwrite(
  bout,
  file = "/Users/lovell/Desktop/manuscripts/genespace_2022/siData/siData4_grassSynBlkCoords.txt",
  sep = "\t", quote = F)

maizePg <- convert_pgdb2wide(
  "/Users/lovell/Desktop/manuscripts/genespace_2022/GENESPACE_data/results/grasses/results/maize_pangenomeDB.txt.gz")
fwrite(
  maizePg,
  file = "/Users/lovell/Desktop/manuscripts/genespace_2022/siData/siData5_maizeAnchoredPangenomeAnnotation.txt",
  sep = "\t", quote = F)


################################################################################
# -- maize, SI data 6-7
blks <- fread("/Users/lovell/Desktop/manuscripts/genespace_2022/GENESPACE_data/results/maize/results/syntenicBlocks.txt.gz", na.strings = c("", "NA"), showProgress = F)
bout <- with(blks, data.table(
  genome1 = gen1, genome2 = gen2, chr1 = chr1, chr2 = chr2,
  start1 = startBp1, end1 = endBp1, start2 = startBp2, end2 = endBp2,
  orient = orient, nHits = nHits1))
fwrite(
  bout,
  file = "/Users/lovell/Desktop/manuscripts/genespace_2022/siData/siData6_maizeNamSynBlkCoords.txt",
  sep = "\t", quote = F)

maizePg <- convert_pgdb2wide(
  "/Users/lovell/Desktop/manuscripts/genespace_2022/GENESPACE_data/results/maize/results/B73_pangenomeDB.txt.gz")
fwrite(
  maizePg,
  file = "/Users/lovell/Desktop/manuscripts/genespace_2022/siData/siData7_maizeAnchoredPangenomeAnnotation.txt",
  sep = "\t", quote = F)

# -- see maize qtl script for generation of SI data 8

################################################################################
# -- Rho duplicate SI data 9-10
blks <- fread("/Users/lovell/Desktop/manuscripts/genespace_2022/GENESPACE_data/results/rho/results/syntenicBlocks.txt.gz", na.strings = c("", "NA"), showProgress = F)
bout <- with(blks, data.table(
  genome1 = gen1, genome2 = gen2, chr1 = chr1, chr2 = chr2,
  start1 = startBp1, end1 = endBp1, start2 = startBp2, end2 = endBp2,
  orient = orient, nHits = nHits1))
fwrite(
  bout,
  file = "/Users/lovell/Desktop/manuscripts/genespace_2022/siData/siData9_rhoGrassSynBlkCoords.txt",
  sep = "\t", quote = F)

hits <- fread("/Users/lovell/Desktop/manuscripts/genespace_2022/GENESPACE_data/results/rho/results/Sviridis_Phallii_synHits.txt.gz", na.strings = c("", "NA"))
hitsog <- subset(hits, isOg & isAnchor & isRep1 & isRep2)
load("/Users/lovell/Desktop/manuscripts/genespace_2022/analysis/Sviridis_Phallii_align.rda")
hitsog[,`:=`(pid1 = sapply(alignList, pid, type = "PID1"),
           pid2 = sapply(alignList, pid, type = "PID2"))]
hitsog[,isRho := grepl("sec", blkID)]
hitsog[,isOverRetain := grepl("prim", blkID) & chr1 == "Chr_08" & chr2 == "Chr03"]
hitsog[,isDefOg := gsub("_", "", chr1) == chr2]
hitsog[,blkcat := ifelse(isRho, "rho", ifelse(isDefOg, "orth", ifelse(isOverRetain, "overr", "ambig")))]

gff <- fread("/Users/lovell/Desktop/manuscripts/genespace_2022/GENESPACE_data/results/rho/results/gffWithOgs.txt.gz", na.strings = c("", "NA"), showProgress = F)
iv <- gff$id; names(iv) <- gff$ofID
hitsog[,`:=`(id1 = iv[ofID1], id2 = iv[ofID2])]
out <- with(hitsog, data.table(
  id1 = id1, id2 = id2, genome1 = gen1, genome2 = gen2, chr1 = chr1, chr2 = chr2,
  start1 = start1, end1 = end1, ord1 = ord1,
  start2 = start2, end2 = end2, ord2 = ord2,
  blkID = blkID, regID = regID, pid1 = pid1, pid2 = pid2, blockType = blkcat))

fwrite(
  out,
  file = "/Users/lovell/Desktop/manuscripts/genespace_2022/siData/siData10_SviridisPhalliiSynAnchorOgPids.txt",
  sep = "\t", quote = F)

