
data <- read.table("SWATH/Cell_Study_SRL_1&3_Re_Analysis_for_ET_180409.csv", row.names=1, header=TRUE, sep=',') # 5854 x 111
stopifnot(data[1,1] == 8677605.061)

cols <- colnames(data) # Sample names for the proteins
cols <- gsub("..sample.1..", "", cols)
cols <- gsub(".wiff", "", cols)
cols <- gsub("SM_PC5", "PC5_SM", cols)
cols <- gsub("X1", "1", cols)
cols <- gsub("Sample", "sample", cols)
cols <- gsub("171230_", "", cols)
cols <- gsub("180107_", "", cols)
cols <- gsub("180110_", "", cols)
cols <- gsub("180111_", "", cols)
cols <- gsub("180112_", "", cols)
cols <- gsub("171217_", "", cols)
cols <- gsub("Study",  "study",  cols)
cols <- gsub("Sample", "sample", cols)
cols <- gsub("sample41", "sample_41", cols)
cols <- gsub("sample43", "sample_43", cols)
cols <- gsub("sample49", "sample_49", cols)
cols <- gsub("sample50", "sample_50", cols)
cols <- gsub("sample51", "sample_51", cols)
cols <- gsub("sample52", "sample_52", cols)
cols <- gsub("sample53", "sample_53", cols)
cols <- gsub("sample54", "sample_54", cols)
cols <- gsub("sample56", "sample_56", cols)
cols <- gsub("sample57", "sample_57", cols)
cols <- gsub("sample58", "sample_58", cols)
cols <- gsub("sample72", "sample_72", cols)
cols <- gsub("sample76", "sample_76", cols)
cols <- gsub("sample80", "sample_80", cols)
cols <- gsub("sample81", "sample_81", cols)
cols <- gsub("sample88", "sample_88", cols)
cols <- gsub("sample90", "sample_90", cols)
cols <- gsub("Cell.study", "Cell_study", cols)
cols <- gsub("Cell study", "Cell_study", cols)
cols <- unlist(lapply(strsplit(as.character(cols), "..", fixed=TRUE), '[[', 2))
colnames(data) <- cols

info <- read.table("SWATH/SWATHDetails.csv", row.names=1, header=TRUE, sep=',')
info$ID <- row.names(info)
info <- info[, c("ID", "Sample", "SWATH.Processing.Date", "Acq_Type", "SWATH.Processing.Instrument", "SWATH.File.Name", "IDA.run.Date", "MS_used", "Operator", "MS_Method")]
info <- info[with(info, order(Sample, SWATH.File.Name)),]
stopifnot(nrow(info) == length(cols))

info$SWATH.File.Name <- gsub("Study",  "study",  info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("Sample", "sample", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("SM_PC5", "PC5_SM", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("171230_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("180107_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("180110_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("180111_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("180112_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("171217_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub(" Cell", "Cell", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub(" SWATH", "SWATH", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample58", "sample_58", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("171230_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("180107_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("180110_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("180111_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("180112_", "", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("SM_PC5", "PC5_SM", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample41", "sample_41", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample43", "sample_43", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample49", "sample_49", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample50", "sample_50", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample51", "sample_51", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample52", "sample_52", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample53", "sample_53", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample54", "sample_54", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample56", "sample_56", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample57", "sample_57", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample58", "sample_58", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample72", "sample_72", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample76", "sample_76", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample80", "sample_80", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample81", "sample_81", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample88", "sample_88", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("sample90", "sample_90", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("Cell.study", "Cell_study", info$SWATH.File.Name)
info$SWATH.File.Name <- gsub("Cell study", "Cell_study", info$SWATH.File.Name)

# Remove "IIICF/a2_A" (IMPORTANT)
data <- data[,c(-1)]; stopifnot(!("171230_SM_PC5_S_Cell_Study_Sample_256" %in% colnames(data)))

map <- merge(data.frame(Cols=cols), info, by.x="Cols", by.y="SWATH.File.Name")
map <- map[with(map, order(Sample)),]
write.table(map$Sample, file="SWATH/sample.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

map$NewSample <- read.table("SWATH/newSample.txt", header=FALSE, sep="\n")$V1
system("python3 SWATH/SWATH_CLEANUP.py")

# Substitue the intensity table with the samples
colnames(data) <- as.character(map$NewSample[match(colnames(data), map$Cols)])

nInfo <- merge(info, map, by.x="ID", by.y="ID")
nInfo <- data.frame(ID=nInfo$ID,
                    RunDate=nInfo$IDA.run.Date.x,
                    Sample=nInfo$NewSample,
                    ProcessDate=nInfo$SWATH.Processing.Date.x,
                    Type=nInfo$Acq_Type.x,
                    Instrument=nInfo$SWATH.Processing.Instrument.x,
                    MSUsed=nInfo$MS_used.x,
                    MSMethod=nInfo$MS_Method.x,
                    Operator=nInfo$Operator.x)

noRep <- function(x) { gsub('_r1', '', gsub('_r2', '', gsub('_r3', '', gsub('_r4', '', x)))) }
nInfo$SampleNR <- noRep(nInfo$Sample)
nInfo <- nInfo[with(nInfo, order(SampleNR)),]

mortals <- c("JFCF_6", "GM02063", "IIICF_E6E7_C4_pre", "IVG_BF_LXSN_pre", "LFS_05F_24_pre", "MeT_4A_pre", "WI38", "IIICF_P7", "IIICF_P9")
nInfo$Mortality <- "Immortal"
nInfo[nInfo$SampleNR %in% mortals,]$Mortality <- "Mortal"

nInfo$Cell <- sapply(strsplit(nInfo$SampleNR, "_"), `[`, 1)
nInfo[nInfo$Cell == "GM02063" | nInfo$Cell == "GM847",]$Cell <- "GM"
nInfo[nInfo$Cell == "VA13" | nInfo$Cell == "WI38",]$Cell <- "VAWI"

# Remove the sample recommended by Erdahl
nInfo <- nInfo[nInfo$ID != 256,]

write.table(data,  file="SWATH/data2.tsv", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(nInfo, file="SWATH/nInfo.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
