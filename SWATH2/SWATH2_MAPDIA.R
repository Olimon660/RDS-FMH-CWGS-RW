
trk   <- read.table("SWATH2/SWATH2_track.tsv", header=1, sep='\t')
data  <- read.table("SWATH2/SWATH2_data.tsv", header=1, sep='\t')
tests <- read.table("SWATH2/SWATH2_tests.tsv", header=1, sep='\t')

write.table(getTests(), file="SWATH2/SWATH2_tests.tsv", row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE)
