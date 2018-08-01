
data <- read.table("/Users/twong/Desktop/SWATH/Cell_Study_SRL_1&3_Re_Analysis_for_ET_180409.csv", row.names=1, header=TRUE, sep=',') # 5854 x 111
cols <- colnames(data) # Sample names for the proteins

dets <- read.table("/Users/twong/Desktop/SWATH/SWATHDetails.csv", row.names=1, header=TRUE, sep=',')
dets <- dets[, c("Sample", "SWATH.File.Name", "SWATH.file.Location")]
stopifnot(nrow(det) == length(cols))

