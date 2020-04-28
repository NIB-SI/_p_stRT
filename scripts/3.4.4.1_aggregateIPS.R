# R script

if (!require(data.table)) install.packages("data.table")
library(data.table)

# assuming running script from ./scripts subdirectory

pathTo = "../genotype_IPS" # change accordingly 

# if okay and okalt sets were sun separately
DT1 = read.delim(file= paste0(pathTo, "/all.okay.aa.tsv"), 
                 header = FALSE, sep = "\t", quote = NULL, 
                 dec = ".", stringsAsFactors = FALSE, 
                 na.strings = "NA", fill = TRUE)
DT2 = read.delim(file= paste0(pathTo, "/all.okalt.aa.tsv"), 
                 header = FALSE, sep = "\t", quote = NULL, 
                 dec = ".", stringsAsFactors = FALSE, 
                 na.strings = "NA", fill = TRUE)

# merge results for okay and okalt
DT = rbind(DT1, DT2)

# fill empty cells with -
DT$V12[DT$V12 == ""] = "-"
DT$V13[DT$V13 == ""] = "-"
DT$V4[DT$V4 == ""] = "-"
DT$V5[DT$V5 == ""] = "-"

# in the case of zero hits, set e-value to a high number
DT$V9[DT$V9 == '-'] = 1.0e+10
DT$V9 = as.numeric(DT$V9)

# keep only resulty with e-value equal or lower that e.g. 1.0e-0 or 1.0e-5
DT3 = DT[DT$V9 <=  1.0e-0,] # 1.0e-5,]
setDT(DT3)

# select which data bases to ignore, if any; if not - skip the next three lines
ignoreDB = c('Coils', 'Gene3D', 'MobiDBLite', 'SUPERFAMILY')
ind = which(DT3$V4 %in% ignoreDB)
DT3 = DT3[-ind,]

# aggregate results
DT4 = DT3[ , .(X = paste(V4, collapse="; ")), by = V1]
DT5 = DT3[ , .(Y = paste(V5, collapse="; ")), by = V1]
DT6 = DT3[ , .(Q = paste(V12, collapse="; ")), by = V1]
DT7 = DT3[ , .(P = paste(V13, collapse="; ")), by = V1]
DT4$X = sapply(1:nrow(DT4), function(x) paste(sort(unique(trimws(unlist(
  strsplit(DT4$X[x],split="; ",fixed=TRUE))))), collapse = "; "))
DT5$Y = sapply(1:nrow(DT5), function(x) paste(sort(unique(trimws(unlist(
  strsplit(DT5$Y[x],split="; ",fixed=TRUE))))), collapse = "; "))
DT6$Q = sapply(1:nrow(DT6), function(x) paste(sort(unique(trimws(unlist(
  strsplit(DT6$Q[x],split="; ",fixed=TRUE))))), collapse = "; "))
DT7$P = sapply(1:nrow(DT7), function(x) paste(sort(unique(trimws(unlist(
  strsplit(DT7$P[x],split="; ",fixed=TRUE))))), collapse = "; "))

DT4$X = sapply(1:nrow(DT4), function(x) gsub("-; ", "", DT4$X[x]))
DT5$Y = sapply(1:nrow(DT5), function(x) gsub("-; ", "", DT5$Y[x]))
DT6$Q = sapply(1:nrow(DT6), function(x) gsub("-; ", "", DT6$Q[x]))
DT7$P = sapply(1:nrow(DT7), function(x) gsub("-; ", "", DT7$P[x]))
total <- merge(DT4,DT5, by="V1")
total <- merge(total,DT6, by="V1")
total <- merge(total,DT7, by="V1")


# see: https://github.com/ebi-pf-team/interproscan/wiki/InterProScan5OutputFormats to define column names
colnames(total) = c("polypeptide_ID", "Analysis", "Signature_Accession", 
                    "IPR_annotations_accession", "IPR_annotations_description")

pathToOutputFolder = "../output" # change accordingly 
# write aggregated results
write.table(total, file = paste0(pathToOutputFolder, 
                                 "/genotype_IPS_filtered_aggregated_filtered.tsv"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", 
            dec = ".", row.names = FALSE, col.names = TRUE)
sessionInfo()

