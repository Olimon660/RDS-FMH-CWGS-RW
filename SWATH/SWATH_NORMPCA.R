
data <- read.table("/Users/twong/Sources/RDS-FMH-CWGS-RW/SWATH/data2.txt", row.names=1, header=TRUE, sep=' ')
print(dim(data))

dets <- read.table("/Users/twong/Sources/RDS-FMH-CWGS-RW/SWATH/newSWATHDetails.tsv", sep='\t')
dets <- dets[with(dets, order(Sample, SWATH.File.Name)),]


