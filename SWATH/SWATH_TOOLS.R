
.remove <- function()
{
    c("Observed.RT.", ".wiff", " ")
}

computerFriendly <- function(x)
{
    for (i in .remove())
        x <- gsub(i, "", x)
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
    track$FileName <- computerFriendly(track$FileName)
    track
}
