# Script to generate Table1 - comparison of run types for cotton
setwd("/Users/lovell/Desktop/manuscripts/genespace_2022/data")
library(GENESPACE)
library(ggplot2)
load("../gsParamList.rda", verbose = T)
gpar <- gsParamList$cotton
gspl <- gsParamList$cottonSplit
rm(list = "gsParamList")

# -- Read in the gffs
gtet <- fread("../cotton_tetraploid/results/gffWithOgs.txt.gz",
              na.strings = c("", "NA"))
gspl <- fread("../cotton_split/results/gffWithOgs.txt.gz",
              na.strings = c("", "NA"))
gsplout <- fread("../cotton_split_outgroup/results/gffWithOgs.txt.gz",
              na.strings = c("", "NA"))
# -- combine (long format)
gtet[,run := "tet"]
gspl[,run := "spl"]
gff <- rbind(gtet, gspl)

# -- subset to non-scaffold and array rep genes
gff <- subset(gff, !grepl("scaff", chr) & isArrayRep)

# -- reformat chr IDs, splitting subgenome and chr numbers
gff[,`:=`(chrn = as.numeric(substr(chr, 2, 3)), subg = substr(chr,1,1))]

# -- subset to non chr2-5
gff <- subset(gff, !chrn %in% 2:5)

# -- strip off A/D genome classes in split run
gff[,genome := gsub("A|D", "", genome)]

# -- reformat to long
m <- melt(
  gff,
  id.vars = c("ofID", "arrayID", "run", "chrn", "subg", "genome"),
  measure.vars = c("globOG", "synOG", "og"),
  variable.name = "runType", value.name = "ogID")

# -- count chrs, subgen, genomes, geneIDs for each run, OG type and OGID
cnts <- m[,list(nChrns = uniqueN(chrn),
                nSubg = uniqueN(subg),
                nIDs = uniqueN(ofID),
                nGenome = uniqueN(genome)),
          by = c("ogID", "runType","run")]

# -- classify as single copy or not
cnts[,is1x := nChrns == 1 & nSubg == 2, nIDs == 2]

# -- count unique OGs that are single copy for each runType and OG
cnt <- cnts[,list(uniqueN(ogID[is1x])),
            by = c("runType","run")]

# -- reformat wide
out <- dcast(cnt, runType ~ run, value.var = "V1")

# -- calculate present better (spl over tetraploid)
out[,percBetter := round(100*((spl - tet)/spl),1)]

# -- print result
knitr::kable(out[,c("runType","tet", "spl","percBetter")])

#   |runType |   tet|   spl| percBetter|
#   |:-------|-----:|-----:|----------:|
#   |globOG  | 15383| 18304|       16.0|
#   |synOG   | 16209| 21530|       24.7|
#   |og      | 22296| 21880|       -1.9|


h <- fread(
  file.path(gpar$paths$results, "Gbarbadense_Gbarbadense_synHits.txt.gz"),
  na.strings = c("NA", ""), showProgress = F)
gff <- fread("../cotton_tetraploid/results/gffWithOgs.txt.gz",
              na.strings = c("", "NA"))

# -- subset to non-scaffold and get chr coords
gff <- subset(gff, !grepl("scaff", chr) & genome == "Gbarbadense")

# -- reformat chr IDs, splitting subgenome and chr numbers
gff[,`:=`(chrn = as.numeric(substr(chr, 2, 3)), subg = substr(chr,1,1))]

# -- subset to non chr2-5
gff <- subset(gff, !chrn %in% 2:5)

# -- get synteny parameters for just Pima cotton
pimaSyn <- subset(gpar$params$synteny, genome1 == genome2 & genome1 == "Gbarbadense")
idsA <- subset(gff, subg == "A")$ofID
idsD <- subset(gff, subg == "D")$ofID
pimaHits <- subset(h, ofID1 %in% idsA & ofID2 %in% idsD)
pimaHits <- subset(pimaHits, substr(chr1, 2, 3) == substr(chr2, 2, 3))
pimaHits[,u := paste(ofID1, ofID2)]

# -- calculate block coordinates from raw hits
tmp <- run_mcscanx(
  gsParam = gpar,
  hits = pimaHits,
  blkSize = pimaSyn$blkSize,
  nGaps = pimaSyn$nGaps,
  path2mcscanx = gpar$paths$mcscanxCall)
rawMcs <- subset(pimaHits, u %in% names(tmp))
rawMcs[,blkID := tmp[u]]
rawMcs <- calc_blkCoords(rawMcs)

# -- calculate block coordinates from collinear array reps
tmp <- run_mcscanx(
  gsParam = gpar,
  hits = subset(pimaHits, isRep1 & isRep2),
  blkSize = pimaSyn$blkSize,
  nGaps = pimaSyn$nGaps,
  path2mcscanx = gpar$paths$mcscanxCall)
repMcs <- subset(pimaHits, u %in% names(tmp))
repMcs[,blkID := tmp[u]]
repMcs <- calc_blkCoords(repMcs)

# -- calculate block coordinates from only hits in the same orthogroup
tmp <- run_mcscanx(
  gsParam = gpar,
  hits = subset(pimaHits, isOg),
  blkSize = pimaSyn$blkSize,
  nGaps = pimaSyn$nGaps,
  path2mcscanx = gpar$paths$mcscanxCall)
ogMcs <- subset(pimaHits, u %in% names(tmp))
ogMcs[,blkID := tmp[u]]
ogMcs <- calc_blkCoords(ogMcs)

# -- pull corresponding syntenic block coords from genespace run
gsBlks <- subset(pimaHits, !is.na(blkID) & isAnchor)
gsBlks <- calc_blkCoords(gsBlks)

# -- reformat blks into single window
gsb <- rbind(
  with(rawMcs, data.table(
    chr = chr1, start = startBp1,
    end = endBp1, blkID = paste0(blkID, "_A"), type = "mcs_raw")),
  with(rawMcs, data.table(
    chr = chr2, start = minBp2,
    end = maxBp2, blkID = paste0(blkID, "_D"), type = "mcs_raw")),
  with(repMcs, data.table(
    chr = chr1, start = startBp1,
    end = endBp1, blkID = paste0(blkID, "_A"), type = "mcs_rep")),
  with(repMcs, data.table(
    chr = chr2, start = minBp2,
    end = maxBp2, blkID = paste0(blkID, "_D"), type = "mcs_rep")),
  with(ogMcs, data.table(
    chr = chr1, start = startBp1,
    end = endBp1, blkID = paste0(blkID, "_A"), type = "mcs_og")),
  with(ogMcs, data.table(
    chr = chr2, start = minBp2,
    end = maxBp2, blkID = paste0(blkID, "_D"), type = "mcs_og")),
  with(gsBlks, data.table(
    chr = chr1, start = startBp1,
    end = endBp1, blkID = paste0(blkID, "_A"), type = "gs")),
  with(gsBlks, data.table(
    chr = chr2, start = minBp2,
    end = maxBp2, blkID = paste0(blkID, "_D"), type = "gs")))
setkey(gsb, chr, start, end)

# -- make 10kb grid for genic region of chrs
tmp <- subset(gff, ofID %in% c(pimaHits$ofID1, pimaHits$ofID2))
chrCoords <- tmp[,list(start = min(start), end = max(end)), by = "chr"]
wind <- 1e4
grd <- chrCoords[,list(
  start = seq(from = start + (wind/2), to = end - (wind/2), by = wind)),
  by = "chr"]
grd[,`:=`(end = start, grdID = 1:.N)]
setkey(grd, chr, start, end)

# -- calculate overlaps
ovlps <- foverlaps(grd, gsb)

# -- count overlaps by grid
cnts <- ovlps[,list(n = uniqueN(blkID[!is.na(blkID)])), by = c("type", "grdID")]
cnts <- dcast(cnts, grdID ~ type, value.var = "n")
cnts[,`NA` := NULL]
cnts[is.na(cnts)] <- 0
cnts <- melt(cnts, id.vars = "grdID", value.name = "n", variable.name = "runType")
cnts[,cat := ifelse(n > 1, "2+", as.character(n))]
sumCnts <- cnts[,list(n = .N), by = c("runType", "cat")]
sumCnts[,perc := round((n/sum(n)) * 100, 1), by = "runType"]

# -- reformat wide
sumCnts[,runType := factor(runType, levels = c("mcs_raw", "mcs_rep", "mcs_og", "gs"))]
out <- dcast(sumCnts, runType ~ cat, value.var = "perc")

# -- print result
knitr::kable(out)

#   |runType |   0|    1|   2+|
#   |:-------|---:|----:|----:|
#   |mcs_raw | 6.5| 79.5| 14.0|
#   |mcs_rep | 6.1| 83.1| 10.7|
#   |mcs_og  | 6.1| 91.3|  2.6|
#   |gs      | 5.6| 93.7|  0.6|
