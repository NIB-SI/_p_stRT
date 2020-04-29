#zagor
# IDs_withFiltering

rm()
gc()

if(!require(data.table)) install.packages("data.table")
library(data.table)
if(!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)
if(!require(grr)) install.packages("grr")
library(grr)

'%ni%' = Negate('%in%')

# assuming you are in scripts subdirectory
# getwd()




# genotype

## IDs, evg headers and seq lenghts

tr = read.table(file = "../intermediate/genotype/genotype.tr_lengths.tsv", 
               header = TRUE, 
               sep = "\t", 
               quote = NULL,
               dec = ".", 
               stringsAsFactors = FALSE,
               comment.char = "@")

cds = read.table(file = "../intermediate/genotype/genotype.cds_lengths.tsv", 
               header = TRUE, 
               sep = "\t", 
               quote = NULL,
               dec = ".", 
               stringsAsFactors = FALSE,
               comment.char = "@")

colnames(tr)[1] = "name"
colnames(cds)[1] = "name"

head(tr)
head(cds)

tr = tr[, -c(2,3)]
cds = cds[, -c(2,3)]

tmp = strsplit(tr$name, " ")
tmp = data.table::transpose(tmp)
length(tmp)
ID = tmp[[1]]
evgHead1 = tmp[[2]]
evgHead2 = tmp[[3]]
## gsub("[;\"]"
evgHead2 = sapply(1:length(evgHead2), function(x) gsub("aalen=", "", evgHead2[x]))
evgHead2 = sapply(1:length(evgHead2), function(x) gsub(";", "", evgHead2[x]))
temp = strsplit(evgHead2, ",")
temp = data.table::transpose(temp)
evgAAlen = temp[[1]]
evgAAperc = temp[[2]]
evgAAcomplete = temp[[3]]

evgHead1 = sapply(1:length(evgHead1), function(x) gsub("evgclass=", "", evgHead1[x]))
evgHead1 = sapply(1:length(evgHead1), function(x) gsub(";", "", evgHead1[x]))
temp = strsplit(evgHead1, ",")
temp = data.table::transpose(temp)
class = temp[[1]]
pass = temp[[2]]
match = temp[[3]]
pct = temp[[4]]

tr.2 = cbind(ID, class, pass, match, pct, evgAAlen, evgAAperc, evgAAcomplete, tr$length)
colnames(tr.2)[dim(tr.2)[2]] = "length"

tmp = strsplit(cds$name, " ")
tmp = data.table::transpose(tmp)
length(tmp)
ID = tmp[[1]]
evgHead1 = tmp[[8]]
evgHead2 = tmp[[3]]

evgHead2 = sapply(1:length(evgHead2), function(x) gsub("aalen=", "", evgHead2[x]))
evgHead2 = sapply(1:length(evgHead2), function(x) gsub(";", "", evgHead2[x]))
temp = strsplit(evgHead2, ",")
temp = data.table::transpose(temp)
evgAAlen = temp[[1]]
evgAAperc = temp[[2]]
evgAAcomplete = temp[[3]]

evgHead1 = sapply(1:length(evgHead1), function(x) gsub("evgclass=", "", evgHead1[x]))
evgHead1 = sapply(1:length(evgHead1), function(x) gsub(";", "", evgHead1[x]))
temp = strsplit(evgHead1, ",")
temp = data.table::transpose(temp)
class = temp[[1]]
pass = temp[[2]]
match = temp[[3]]
pct = temp[[4]]

clen = gsub("clen=", "", tmp[[4]])
strand = gsub("strand=", "", tmp[[5]])
offs = gsub("offs=", "", tmp[[6]])
clen = gsub(";", "", clen)
strand = gsub(";", "", strand)
offs = gsub(";", "", offs)

# dim(cbind(clen, strand, offs))

cds.2 = cbind(ID, class, pass, match, pct, evgAAlen, evgAAperc, evgAAcomplete, clen, strand, offs, cds$length)
colnames(cds.2)[dim(cds.2)[2]] = "length"


colnames(tr.2)[2:dim(tr.2)[2]] = paste0("tr.", colnames(tr.2)[2:dim(tr.2)[2]])
colnames(cds.2)[2:dim(cds.2)[2]] = paste0("cds.", colnames(cds.2)[2:dim(cds.2)[2]])

tr.3 = merge(tr.2, cds.2, by="ID", all.x =TRUE)

tr.3$tr.length = as.numeric(levels(tr.3$tr.length))[tr.3$tr.length]
tr.3$cds.length = as.numeric(levels(tr.3$cds.length))[tr.3$cds.length]




## lenghts

hist(as.numeric(tr.3$tr.length), 
     breaks = seq(0, max(as.numeric(tr.3$tr.length), na.rm = TRUE) + 50, 50),
     main = "genotype evigene transcripts length", xlab = "transcript length",
     xlim=c(0, max(300,max(as.numeric(tr.3$tr.length), na.rm = TRUE))),
     col = "grey90")
hist(as.numeric(tr.3$cds.length), 
     breaks = seq(0, max(as.numeric(tr.3$cds.length), na.rm = TRUE) + 50, 50),
     main = "genotype evigene coding sequences length", xlab = "CDS length",
     xlim=c(0, max(300,max(as.numeric(tr.3$cds.length), na.rm = TRUE))),
     col = "grey45", las=2)
abline(v=180, col = "red")
abline(v=200, col = "blue")
abline(v=300, col = "green")
axis(1, at=c(180), labels=c(180), las=2, col = c("red"), cex.axis=0.75, col.axis = c("red"))
axis(1, at=c(200), labels=c(200), las=2, col = c("blue"), cex.axis=0.75, col.axis = c("blue"))
axis(1, at=c(300), labels=c(300), las=2, col = c("green"), cex.axis=0.75, col.axis = c("green"))





## Match Annot

MA = read.table(file = "../intermediate/genotype/genotype.tr.okay_ref_STARlong.matchAnnot.parsed.txt", 
               header = FALSE, 
               sep = " ", 
               quote = NULL,
               dec = ".", 
               stringsAsFactors = FALSE,
               na.strings = "",
               comment.char = "@")

# dim(MA)
head(MA)
colnames(MA) = c("ID", "DMgeneID", "DMtrID", "exon_match", "match_score")
MA = MA[order(MA$ID),] 
duplRow = which(duplicated(MA) | duplicated(MA[nrow(MA):1, ])[nrow(MA):1])
MA[duplRow[1:10],]

MA = MA[!duplicated(MA),]

hist(MA$match_score, main = 'All Match Annot scores', xlab = "MA score")


data <- data.frame(cbind(seq(1,dim(MA)[1],1),MA))
colnames(data)[1] = "enumerate"
setDT(data)
DT5 = data[ , .(match.score = paste(match_score,
                         collapse=",")), by = ID]

no.gene_score.0 = 
    unlist(sapply(1:length(DT5$match.score), function(i) {
      aa = as.numeric(unlist(strsplit(DT5$match.score[i], ",")))
      if(all(is.na(aa)) | max(aa, na.rm=TRUE) == 0) return(i)
    }))
ind = unlist(matches((unlist(DT5[no.gene_score.0 ,1])), MA[,1],  list =TRUE))
no.gene_score.0 = MA[ind ,]
unique(no.gene_score.0[,5])

MA.2 = MA[-which(MA$match_score == 0),]
MA = MA.2[-which(is.na(MA.2$match_score)),]

# only ones with max score
# uniqueID = unique(MA$ID)
keep = NULL

data <- data.frame(cbind(seq(1,dim(MA)[1],1),MA))
colnames(data)[1] = "enumerate"

dd = data %>% group_by(ID) %>% summarize(maxID = paste(enumerate[which(match_score == max(match_score))],  
                                                       collapse=","))



keep = as.numeric(unlist(strsplit(dd$maxID, ",")))
dim(MA)[1] - length(keep)
MA = MA[keep,]
hist(MA$match_score, main = 'Selected Match Annot scores', xlab = "MA score")

setDT(MA)
setDT(no.gene_score.0)
# dim(MA)
DT1 = MA[ , .(DM.geneID = paste(DMgeneID,
                         collapse="; ")), by = ID]
DT2 = MA[ , .(DM.trID = paste(DMtrID,
                         collapse="; ")), by = ID]
DT3 = MA[ , .(exon.match = paste(exon_match,
                         collapse="; ")), by = ID]
DT4 = MA[ , .(match.score = paste(match_score,
                         collapse="; ")), by = ID]
DT1.1 = no.gene_score.0[ , .(DM.geneID = paste(DMgeneID,
                         collapse="; ")), by = ID]
DT2.1 = no.gene_score.0[ , .(DM.trID = paste(DMtrID,
                         collapse="; ")), by = ID]
DT3.1 = no.gene_score.0[ , .(exon.match = paste(exon_match,
                         collapse="; ")), by = ID]
DT4.1 = no.gene_score.0[ , .(match.score = paste(match_score,
                         collapse="; ")), by = ID]

# dim(MA)

myVec = unique(MA[,1])
myVec.1 = unique(no.gene_score.0[,1])
# dim(myVec)
# dim(myVec.1)
# dim(DT1)
# dim(DT1.1)

myVec = merge(myVec, DT1, by="ID", all.x =TRUE)
myVec = merge(myVec, DT2, by="ID", all.x =TRUE)
myVec = merge(myVec, DT3, by="ID", all.x =TRUE)
myVec = merge(myVec, DT4, by="ID", all.x =TRUE)
myVec.1 = merge(myVec.1, DT1.1, by="ID", all.x =TRUE)
myVec.1 = merge(myVec.1, DT2.1, by="ID", all.x =TRUE)
myVec.1 = merge(myVec.1, DT3.1, by="ID", all.x =TRUE)
myVec.1 = merge(myVec.1, DT4.1, by="ID", all.x =TRUE)

# dim(myVec)
# dim(myVec.1)

no.gene_score.0 = merge(tr.3, myVec.1, all.y =TRUE)

# print(dim(merge(tr.3, myVec, by="ID", all.x =TRUE)))

tr.3 = merge(tr.3, myVec, by="ID", all.x =TRUE)
colnames(myVec)[1] = "ID2"






## IPR

IPR = read.table(file = "../intermediate/genotype/genotype_IPS_filtered_aggregated_filtered.tsv", 
               header = TRUE, 
               sep = "\t", 
               quote = NULL,
               dec = ".", 
               stringsAsFactors = FALSE,
               comment.char = "@")

# dim(IPR)
colnames(IPR)[1] = c("ID")

# print(dim(merge(tr.3, IPR, by="ID", all.x =TRUE)))

tr.3 = merge(tr.3, IPR, by="ID", all.x = TRUE)

no.gene_score.0 = merge(no.gene_score.0, IPR, all.x = TRUE)

colnames(tr.3)[which(colnames(tr.3) == "Analysis")] = "IPS_Analysis"
colnames(no.gene_score.0)[which(colnames(no.gene_score.0) == "Analysis")] = "IPS_Analysis"

colnames(tr.3)[which(colnames(tr.3) == "Signature_Accession")] = "IPS_Signature_Accession"
colnames(no.gene_score.0)[which(colnames(no.gene_score.0) == "Signature_Accession")] = "IPS_Signature_Accession"




## vecscreen

vs = read.table(file = "../intermediate/genotype/genotype_vecscreen.tsv", 
               header = TRUE, 
               sep = "\t", 
               quote = NULL,
               dec = ".", 
               stringsAsFactors = FALSE,
               comment.char = "@")

# dim(vs)
colnames(vs)[1] = "ID"
length(unique(vs$ID))

setDT(vs)
# dim(vs)
DT1 = vs[ , .(coverage = paste(paste0(Matching_vector._starting_with_uv., " [",
                                       Lower_end_of_the_alignment_in_the_vector, "-",
                                       Upper_end_of_the_alignment_in_the_vector, "] ",
                                      coverage, "%"),
                         collapse="; ")), by = ID]
DT2 = vs[ , .(vectorEvidence = paste(paste0(The_strength_of_this_vecscreen_match_, ", ",
                                       The_strength_of_the_strongest_vecscreen_match_for_this_query, ", ",
                                       Whether_there_is_any_dangling_part_.called_.Suspect._by_vecscreen._at_either_end_of_the_query, ", ",
                                       the_classification_of_the_match),
                         collapse="; ")), by = ID]
DT3 = vs[ , .(blastnVector = paste(blastn_desc, collapse="; ")), by = ID]


myVec = unique(vs[,1])
# dim(myVec)
# dim(DT1)

myVec = merge(myVec, DT1, by="ID", all.x =TRUE)
myVec = merge(myVec, DT2, by="ID", all.x =TRUE)
myVec = merge(myVec, DT3, by="ID", all.x =TRUE)


# print(dim(merge(tr.3, myVec, by="ID", all.x =TRUE)))

tr.3 = merge(tr.3, myVec, by="ID", all.x =TRUE)
no.gene_score.0 = merge(no.gene_score.0, myVec, all.x = TRUE)
colnames(myVec)[1] = "ID2"


colnames(tr.3)[which(colnames(tr.3) == "coverage")] = "VecScreen_coverage"
colnames(no.gene_score.0)[which(colnames(no.gene_score.0) == "coverage")] = "VecScreen_coverage"



## DIAMOND BLASTX

# DIAMOND top hit files renamed: concatenate integer with leading zero at the beggining of the file name

(myfiles = list.files(path = "../intermediate/genotype/", pattern = "_top1.tsv"))

# watch our for comment.chars # in tables
for (i in myfiles) {
  print(i)
  blast = read.table(file = paste0("../intermediate/genotype/",i), 
               header = TRUE, 
               sep = "\t", 
               quote = NULL,
               dec = ".", 
               stringsAsFactors = FALSE,
               comment.char = "@")
  
  blast = blast[,-1]
  tmp = strsplit(blast[,1], " ")
  tmp = data.table::transpose(tmp)[[1]]
  blast[,1] = tmp
  remove = which(duplicated(blast[,1])) # duplicated
  if (length(remove)) blast = blast[-remove,]
  
  colnames(blast)[1] = "ID"
  blast$Gaps.num = round(blast$Gaps.*blast$Aligned_seq_length/100)
  blast$Query_coverage = (blast$Aligned_seq_length-blast$Gaps.num)/(blast$Query_length/3)
  blast$Target_coverage = (blast$Aligned_seq_length-blast$Gaps.num)/blast$Target_length
  blast$Query_Target_ratio = (blast$Query_length/3)/blast$Target_length
  ind1 = which(blast$Query_coverage >= 0.50)
  ind2 = which(blast$Target_coverage >= 0.50)
  ind = intersect(ind1, ind2)
  blast = blast[ind, ]
  blast$Target_coverage[blast$Target_coverage > 1.0] = 1.0
  blast$Query_coverage[blast$Query_coverage > 1.0] = 1.0
  ind = which(colnames(blast) %in% c("ID", "Target_ID", "Target_description", "Aligned_seq_length", "Target_coverage", "Query_coverage", 
                         "E.value", "Score"))
  blast = blast[, ind]
  
  j = gsub(".out.ENCHformat_top1.tsv", "", i)
  j = gsub("_genotype", "", j)
  colnames(blast)[2:dim(blast)[2]] = paste0(j, "_", colnames(blast)[2:dim(blast)[2]])
  
  tr.3 = merge(tr.3, blast, by="ID", all.x =TRUE)
  no.gene_score.0 = merge(no.gene_score.0, blast, by="ID", all.x = TRUE)
}




### save.image("../output/genotype/genotype_withFiltering.RData")
save(tr.3, no.gene_score.0, file = "../output/genotype/tr3s.RData")
rm(list=ls())
gc()



load("../output/genotype/tr3s.RData")

library(data.table)
library(tidyverse)
library(grr)





## write

i <- sapply(tr.3, is.factor) 
tr.3[i] <- lapply(tr.3[i], as.character) 
tr.3[is.na(tr.3)] = NA
tr.3[tr.3 == ""] = NA

i <- sapply(no.gene_score.0, is.factor) 
no.gene_score.0[i] <- lapply(no.gene_score.0[i], as.character) 
no.gene_score.0[is.na(no.gene_score.0)] = NA
no.gene_score.0[no.gene_score.0 == ""] = NA
# write.table(no.gene_score.0, "../output/genotype/no.gene_score.0_all.tsv", 
#             sep = "\t", row.names = FALSE,
#             quote = FALSE)

hist(as.numeric(no.gene_score.0$cds.length), 
     breaks = seq(0, max(as.numeric(no.gene_score.0$cds.length), na.rm = TRUE) + 50, 50),
     main = "genotype evigene MA score NA/0 CDS sequence length", xlab = "CDS length",
     xlim=c(0, max(300,max(as.numeric(no.gene_score.0$cds.length), na.rm = TRUE))),
     col = "grey45", las=2)
abline(v=180, col = "red")
abline(v=200, col = "blue")
abline(v=300, col = "green")
axis(1, at=c(180), labels=c(180), las=2, col = c("red"), cex.axis=0.75, col.axis = c("red"))
axis(1, at=c(200), labels=c(200), las=2, col = c("blue"), cex.axis=0.75, col.axis = c("blue"))
axis(1, at=c(300), labels=c(300), las=2, col = c("green"), cex.axis=0.75, col.axis = c("green"))

table(as.numeric(no.gene_score.0$cds.length) >= 300)
table(as.numeric(no.gene_score.0$cds.length) >= 200)
table(as.numeric(no.gene_score.0$cds.length) >= 180)





## keep or remove

# colnames(tr.3)

okay = which(!is.na(tr.3$DM.trID))
length(okay)
okay = c(okay, which(!is.na(tr.3$IPS_Analysis)))
length(okay)

ind = grep("Target_ID", colnames(tr.3))

okay = c(okay, 
         which(!is.na(tr.3[,ind[1]])),
         which(!is.na(tr.3[,ind[2]])),
         which(!is.na(tr.3[,ind[3]])),
         which(!is.na(tr.3[,ind[4]])),
         which(!is.na(tr.3[,ind[5]])),
         which(!is.na(tr.3[,ind[6]])))

# which(as.numeric(tr.3$cds.length) >= 180)
# which(as.numeric(tr.3$cds.length) < 180)

okay = sort(unique(okay))
length(okay)

tr.3.annot = tr.3[okay,]
tr.3.unknown = tr.3[-okay,]

tr.3.annot$status="keep"
tr.3.unknown$status="discard"

hist(as.numeric(tr.3.annot$tr.length), 
     breaks = seq(0, max(as.numeric(tr.3.annot$tr.length), na.rm = TRUE) + 50, 50),
      main = "OK genotype evigene transcript length", xlab = "transcript length",
     xlim=c(0, max(300,max(as.numeric(tr.3.annot$tr.length), na.rm = TRUE))))
hist(as.numeric(tr.3.annot$cds.length), 
     breaks = seq(0, max(as.numeric(tr.3.annot$cds.length), na.rm = TRUE) + 50, 50),
     main = "OK genotype evigene coding sequence length", xlab = "CDS length",
     xlim=c(0, max(300,max(as.numeric(tr.3.annot$cds.length), na.rm = TRUE))))
abline(v=180, col = "red")
abline(v=200, col = "blue")
abline(v=300, col = "green")
table(as.numeric(tr.3.annot$cds.length) >= 300)
table(as.numeric(tr.3.annot$cds.length) >= 200)
table(as.numeric(tr.3.annot$cds.length) >= 180)

hist(as.numeric(tr.3.unknown$tr.length), 
     breaks = seq(0, max(as.numeric(tr.3.unknown$tr.length), na.rm = TRUE) + 50, 50),
      main = "notOK genotype evigene transcript with secondary ORF length", xlab = "transcript length",
     xlim=c(0, max(300,max(as.numeric(tr.3.unknown$tr.length), na.rm = TRUE))))
hist(as.numeric(tr.3.unknown$cds.length), 
     breaks = seq(0, max(as.numeric(tr.3.unknown$cds.length), na.rm = TRUE) + 50, 50),
     main = "notOK genotype evigene secondary ORF coding sequence length", xlab = "CDS length",
     xlim=c(0, max(300,max(as.numeric(tr.3.unknown$cds.length), na.rm = TRUE))),     
     col = "grey45", las=2)
abline(v=180, col = "red")
abline(v=200, col = "blue")
abline(v=300, col = "green")
axis(1, at=c(180), labels=c(180), las=2, col = c("red"), cex.axis=0.75, col.axis = c("red"))
axis(1, at=c(200), labels=c(200), las=2, col = c("blue"), cex.axis=0.75, col.axis = c("blue"))
axis(1, at=c(300), labels=c(300), las=2, col = c("green"), cex.axis=0.75, col.axis = c("green"))

table(as.numeric(tr.3.unknown$cds.length) >= 300)
table(as.numeric(tr.3.unknown$cds.length) >= 200)
table(as.numeric(tr.3.unknown$cds.length) >= 180)



tr.3.annot[is.na(tr.3.annot)] = "-"
tr.3.unknown[is.na(tr.3.unknown)] = "-"

# write.table(tr.3.annot, file = "../output/genotype/keep/genotype_tr.3_annot.tsv", 
#             append = FALSE, quote = FALSE, sep = "\t",
#             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
#             col.names = TRUE, qmethod = c("escape", "double"),
#             fileEncoding = "")
# 
# write.table(tr.3.unknown, file = "../output/genotype/remove/genotype_tr.3_unknown.tsv", 
#             append = FALSE, quote = FALSE, sep = "\t",
#             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
#             col.names = TRUE, qmethod = c("escape", "double"),
#             fileEncoding = "")

write.table(tr.3.annot$ID, file = "../output/genotype/genotype_keep_IDs.tsv", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(tr.3.unknown$ID, file = "../output/genotype/genotype_drop_IDs.tsv", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")






## combo table

tr.3.unknown$ID2 = tr.3.unknown$ID
tr.3.annot$ID2 = tr.3.annot$ID
combo = rbind(tr.3.unknown, tr.3.annot)
combo = combo[, c(1, dim(combo)[2], seq(2, dim(combo)[2]-1, 1))]
colnames(combo)[1:2]=c("tr.ID", "cds.ID")

write.table(combo, file = "../output/genotype/genotype_tr.cds.tsv", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")



sessionInfo()

rm(list=ls())
gc()
