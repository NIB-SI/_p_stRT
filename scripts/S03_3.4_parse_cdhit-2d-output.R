# zagor
# parsing cdhit-2d assignments

if (!require("data.table")) install.packages("data.table")

# https://stackoverflow.com/questions/13365732/sorting-non-alphanumeric-characters-in-ascii-order-in-r
asciiSort <- function(vec) {
  x <- sapply(vec, 
              function(X) {
                paste0(strtoi(charToRaw(X), base=16), collapse="")
              })
  vec[order(x)]
}


# input file
pathToInput = '../minimal_example_out/cdhit-2d_out' # change to appropriate relative path
pathToOutput = '../output' # change to appropriate relative path
cdhitOut = read.table(file = paste0(pathToInput, '/', "cdhit-2d.clstr"), 
                 header = FALSE, sep = "\n", 
                 quote = NULL,
                 dec = ".", 
                 stringsAsFactors = FALSE, 
                 na.strings = "NA", 
                 comment.char = "#")



cdhitOut = data.table(setDT(cdhitOut))


cdhitOut$V1 = gsub(">Cluster ", "Cluster\t",cdhitOut$V1)
pattern = "\t"
cdhitOut[, c("V1", "V2") := tstrsplit(V1, pattern, fixed=TRUE)]
pattern = "aa, >"
cdhitOut[, c("V2", "V3") := tstrsplit(V2, pattern, fixed=TRUE)]
pattern = "... "
cdhitOut[, c("V3", "V4") := tstrsplit(V3, pattern, fixed=TRUE)]
cdhitOut$V4 = gsub("/", "\t", cdhitOut$V4)

pattern = "\t"
cdhitOut[, c("V4", "V5") := tstrsplit(V4, pattern, fixed=TRUE)]
cdhitOut$V7 = NA
ind = grep("Cluster", cdhitOut$V1)
cdhitOut$V7[ind] = paste0(cdhitOut$V1[ind], " ", cdhitOut$V2[ind])



goodIdx <- !is.na(cdhitOut$V7)
goodVals <- c(NA, cdhitOut$V7[goodIdx])
fillIdx <- cumsum(goodIdx)+1
cdhitOut$V7 = goodVals[fillIdx]

ind = grep("Cluster", cdhitOut$V1)
cdhitOut = cdhitOut[-ind, ]

tmp = table(cdhitOut$V7)
ind = match(cdhitOut$V7, names(tmp))
cdhitOut$V9 = NA
cdhitOut$V9 = tmp[ind]




singletons = cdhitOut[intersect(which(cdhitOut$V9 <= 1), which(cdhitOut$V4 == "*")),]
unique(singletons$V4)
ss = cdhitOut$V7[which(cdhitOut$V9 == 1)]
all (ss %in% singletons$V7)

cdhitOut$V8 = NA
cdhitOut = cdhitOut[order(V7,V4)]
ind = c(grep("\\*", cdhitOut$V4))
cdhitOut$V8[ind] = cdhitOut$V3[ind]


goodIdx <- !is.na(cdhitOut$V8)
any(which(goodIdx) != ind)
goodVals <- c(NA, cdhitOut$V8[goodIdx])
any(is.na(cdhitOut$V8[goodIdx]))
fillIdx <- cumsum(goodIdx) + 1
cdhitOut$V8 = goodVals[fillIdx]

sum (cdhitOut$V7 %in% ss)


representatives = cdhitOut[intersect(which(cdhitOut$V9 > 1), which(cdhitOut$V4 == "*")),]

cdhitOut = cdhitOut[-ind, ]

sum (cdhitOut$V7 %in% ss)

min(cdhitOut$V1)
max(cdhitOut$V1)

colnames(cdhitOut) = c("seqNo", "seqLen", "seqID", "pos", "perc", "cluster", "clu_cnt", "rep")
colnames(representatives) = c("seqNo", "seqLen", "seqID", "pos", "perc", "cluster", "clu_cnt", "rep")
colnames(singletons) = c("seqNo", "seqLen", "seqID", "pos", "perc", "cluster", "clu_cnt")
singletons$rep = singletons$seqID
singletons$wc = NA
singletons$wc[c(grep("PGSC", singletons$seqID), grep("Sotub", singletons$seqID))] = "DM"
singletons$wc[grep('genotype1', singletons$seqID)] = "g1"
singletons$wc[grep('genotype2', singletons$seqID)] = "g2"


cdhitOut$cv.seq = NA
cdhitOut$cv.seq[c(grep("PGSC", cdhitOut$seqID), grep("Sotub", cdhitOut$seqID))] = "DM"
cdhitOut$cv.seq[grep('genotype1', cdhitOut$seqID)] = "g1"
cdhitOut$cv.seq[grep('genotype2', cdhitOut$seqID)] = "g2"


cdhitOut$cv.rep = NA
cdhitOut$cv.rep[c(grep("PGSC", cdhitOut$rep), grep("Sotub", cdhitOut$rep))] = "DM"
cdhitOut$cv.rep[grep('genotype1', cdhitOut$rep)] = "g1"
cdhitOut$cv.rep[grep('genotype2', cdhitOut$rep)] = "g2"



representatives$cv.rep = NA
representatives$cv.seq = NA
representatives$cv.seq = NA
representatives$cv.seq[c(grep("PGSC", representatives$seqID), grep("Sotub", representatives$seqID))] = "DM"
representatives$cv.seq[grep('genotype1', representatives$seqID)] = "g1"
representatives$cv.seq[grep('genotype2', representatives$seqID)] = "g2"


representatives$cv.rep = NA
representatives$cv.rep[c(grep("PGSC", representatives$rep), grep("Sotub", representatives$rep))] = "DM"
representatives$cv.rep[grep('genotype1', representatives$rep)] = "g1"
representatives$cv.rep[grep('genotype2', representatives$rep)] = "g2"



table(cdhitOut$cv.rep)
table(cdhitOut$cv.seq)
table(singletons$wc)

length(cdhitOut$seqID) == length(unique(cdhitOut$seqID))
length(cdhitOut$rep) == length(unique(cdhitOut$rep))
length(singletons$rep) == length(unique(singletons$rep))

length(unique(singletons$rep)) + length(unique(cdhitOut$rep)) + length(unique(cdhitOut$seqID))



a = tapply(cdhitOut$cv.seq, cdhitOut$cluster, FUN=paste)

cdhitOut$combo = sapply(1:length(cdhitOut$cluster), function(x) {
  ind = match(cdhitOut$cluster[x], names(a))
  paste(unique(asciiSort(c(a[[ind]], cdhitOut$cv.rep[x]))), collapse = "|")
})

any(is.na(cdhitOut$combo))

representatives$combo = NA
ind = match(representatives$seqID, cdhitOut$rep)
representatives$combo = cdhitOut$combo[ind]

colnames(singletons)[dim(singletons)[2]] = "cv.seq"
singletons$cv.rep = singletons$cv.seq

singletons$combo = singletons$cv.rep
singletons$status = "singleton"
cdhitOut$status = "alternative"
representatives$status = "representative"


mymerge.tmp.t = rbind(representatives, cdhitOut)
mymerge.tmp.t = rbind(mymerge.tmp.t, singletons)
mymerge = mymerge.tmp.t


write.table(x = mymerge, file = paste0(pathToOutput, '/', "my_cdhit2d_clusters.tsv"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")





