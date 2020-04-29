#zagor
# paralogue clusters from EvidentialGene headers intersected with filtering output

library(stringr)
library(igraph)
library(data.table)

pathToInput = "../input"
pathToOutput = "../output"

data = read.delim(file = paste0(pathToInput, "/IDs.txt"),
                  header = FALSE, sep = "\t", 
                  quote = NULL,
                  dec = ".", 
                  stringsAsFactors = FALSE, 
                  na.strings = "NA", 
                  comment.char = "#",
                  col.names = paste0("ID",seq_len(3)),
                  fill = TRUE)

data[,1] = gsub("utrorf", "", data[,1])
data[,2] = gsub("utrorf", "", data[,2])
data[,3] = gsub("utrorf", "", data[,3])

dataLen = read.delim(file = paste0(pathToInput, "/seq_len.txt"),
                     header = FALSE, sep = "\t", 
                     quote = NULL,
                     dec = ".", 
                     stringsAsFactors = FALSE, 
                     na.strings = "NA", 
                     comment.char = "#",
                     col.names = c("ID", "len"),
                     fill = TRUE)

keepStatus = read.delim(file = paste0(pathToInput, "/genotype_keep_IDs.tsv"),
                     header = FALSE, sep = "\t", 
                     quote = NULL,
                     dec = ".", 
                     stringsAsFactors = FALSE, 
                     na.strings = "NA", 
                     comment.char = "#",
                     col.names = c("ID"),
                     fill = TRUE)


a = data[,1:2]
b = data[,2:3]
colnames(a) = colnames(b) = c("ID1", "ID2")
edges = rbind(a,b)
ind = which(!(edges$ID1 %in% keepStatus$ID) & !(edges$ID2 %in% keepStatus$ID))
edges = edges[-ind,]
edges$ID1[!(edges$ID1 %in% keepStatus$ID)] = ''
edges$ID2[!(edges$ID2 %in% keepStatus$ID)] = ''
edges = edges[-which(edges$ID1== ''),]
edges[which(edges$ID2== ''),]$ID2 = edges[which(edges$ID2== ''),]$ID1
edges = as.matrix(edges)
g = graph_from_edgelist(edges, directed = FALSE)
V(g)
sum(which_multiple(g))
sum(which_loop(g))
g = simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
sum(which_multiple(g))
sum(which_loop(g))
V(g)

comp = components(g, mode = "weak")
table(comp$csize)

clusters = cbind(names(comp$membership), comp$membership)
colnames(clusters) = c("ID", "clu")

clusters = data.table(clusters)
clusters$clu = as.numeric(clusters$clu)
setorder(clusters, clu, ID)


clusters = merge(x = clusters, y = dataLen, by = "ID", all.x=TRUE)
clusters = clusters[!is.na(clusters$len),]
clusters = clusters[order(clusters$clu, clusters$len),]

clusters$status = "alternative"
clusters$alt = clusters$ID
clusters$rep = NA

maxLen = aggregate(len~clu,data=clusters, FUN = max)
ind = match(clusters$clu, maxLen$clu)
clusters$maxLen = maxLen$len[ind]
setorder(clusters, clu, -len)

ind = which(clusters$len == clusters$maxLen)
clusters$status[ind] = "representative"
a = which(duplicated(cbind(clusters$clu, clusters$len)) & (clusters$len == clusters$maxLen))
b = clusters[a,]$clu
clusters$status[a] = "alternative"

clusters$rep = as.character(clusters$rep)
clusters[clusters$status == "representative"]$rep = clusters[clusters$status == "representative"]$alt
clusters[clusters$status != "representative"]$rep = NA
dim(clusters) == dim(clusters[!duplicated(clusters), ])
goodIdx <- !is.na(clusters$rep)
goodVals <- c(NA, clusters$rep[goodIdx])
fillIdx <- cumsum(goodIdx)+1
clusters$rep = goodVals[fillIdx]



# adapt patterns according to your input sequence naming for evigene tr2aacds step, 
#  i.e. assemblerNameGenotypeCount_sequenceCount
clusters$cv.seq[grep("dnRY", clusters$alt)] = "R"
clusters$cv.rep[grep("dnRY", clusters$rep)] = "R"
# or just type genotype abbreviation


n = str_pad(clusters$clu, 6, pad = "0")
clu = paste0("stCusTr_", clusters$cv.rep, "_", n)
clusters$clu = clu

tmp = table(clusters$clu)
ind = match(clusters$clu, names(tmp))
clusters$seqCount = tmp[ind]

clusters.out = clusters[, c(2,9,10,1,4)]

write.table(x = clusters.out, file = paste0(pathToOutput, "/clusters.tsv"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

alt = clusters.out[clusters.out$status == 'alternative',]
rep = clusters.out[clusters.out$status == 'representative',]

write.table(x = alt$ID, 
            file = paste0(pathToOutput, "/alternativesIDs.tsv"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
write.table(x = rep$ID, 
            file = paste0(pathToOutput, "/representativesIDs.tsv"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

alias.rep = cbind(rep$ID, 
              paste0(rep$ID,
                     " [cluID: ",
                     rep$clu,
                     ", cluSize: ",
                     rep$seqCount,
                     "]"
                )
              )
alias.alt = cbind(alt$ID, 
                  paste0(alt$ID,
                         " [cluID: ",
                         alt$clu,
                         ", cluSize: ",
                         alt$seqCount,
                         ", repSeq: ",
                         rep[match(alt$clu, rep$clu),]$ID,
                         "]"
                  )
                )
write.table(x = alias.alt, 
            file = paste0(pathToOutput, "/alias.alt"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
write.table(x = alias.rep, 
            file = paste0(pathToOutput, "/alias.rep"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
