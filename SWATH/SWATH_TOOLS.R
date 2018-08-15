
.remove <- function()
{
    c("Observed.RT.", ".wiff", " ")
}

computerFriendly <- function(x)
{
    for (i in .remove())
        x <- gsub(i, "", x)
    x <- unlist(lapply(strsplit(as.character(x), "..", fixed=TRUE), '[[', 1))    
    x <- gsub("Sample0", "SWATH_Sample_0", x)    
    x <- gsub("Sample1", "SWATH_Sample_1", x)
    x <- gsub("Sample2", "SWATH_Sample_2", x)
    x <- gsub("Sample3", "SWATH_Sample_3", x)
    x <- gsub(".Sample.", "Cell_Study_SWATH_Sample", x)
    x <- gsub("Cell_Study_SWATH_SWATHCell", "Cell", x)
    x <- gsub("Cell_StudyCell_Study", "Cell_Study", x)
    x <- gsub("Cell_studyCell_Study", "Cell_Study", x)
    x <- gsub("CellStudyCell_Study", "Cell_Study", x)
    x <- gsub(".sample_", "Cell_Study_SWATH_Sample", x)
    x
}

track <- function()
{
    track <- read.table("SWATH/Cell Study_Tracking_SM.csv", header=FALSE, sep=',')
    track <- track[,c(1:15)]
    colnames(track) <- c("ID", "Sample", "ProcessDate", "AcqType", "ProcessInst", "FileName", "Location", "Notes", "RunDate", "IDAFile", "IDALocation",
                         "MSUsed", "Operator", "MSMethod", "AcqType")
    track <- track[, c("ID", "Sample", "ProcessDate", "AcqType", "ProcessInst", "FileName", "RunDate", "MSUsed", "Operator", "MSMethod")]
    track <- track[with(track, order(Sample, FileName)),]
    track <- track[track$ID != "",]
    track <- track[track$ID != "Sample No ",]
    track <- track[track$FileName != "",]
    track
}
