
info <- read.table("SWATH2/SWATH2_track.tsv", header=1, sep='\t')
data <- read.table("SWATH2/SWATH2_data.tsv", header=1, sep='\t', check.names=F)
test <- read.table("SWATH2/SWATH2_tests.tsv", header=1, sep='\t')

for (row in 1:nrow(test))
{
    t <- test[row,]
    
    print(t)
    
    
    break
}
