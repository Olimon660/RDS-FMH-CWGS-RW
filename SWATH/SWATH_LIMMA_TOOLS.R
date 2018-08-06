library(limma)

eb.fit <- function(dat, design)
{
    n <- dim(dat)[1]
    fit <- lmFit(dat, design)
    fit.eb <- eBayes(fit)
    logFC <- fit.eb$coefficients[, 2]
    df.r <- fit.eb$df.residual
    df.0 <- rep(fit.eb$df.prior, n)
    s2.0 <- rep(fit.eb$s2.prior, n)
    s2 <- (fit.eb$sigma)^2
    s2.post <- fit.eb$s2.post
    t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2]
    t.mod <- fit.eb$t[, 2]
    p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
    p.mod <- fit.eb$p.value[, 2]
    q.ord <- qvalue(p.ord)$q
    q.mod <- qvalue(p.mod)$q
    results.eb <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
    results.eb[order(results.eb$p.mod), ]
}

getTests <- function()
{
    tests <- read.table("SWATH/contrasts.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
    tests[tests$Mortal == "IIICF",]$Mortal <- "IIICF_P7" # IIICF_P7 is being used as the reference
    tests[tests$Mortal == "JFCF6_P12",]$Mortal <- "JFCF_6"
    tests    
}

getPData <- function()
{
    pData <- data.frame(row.names=colnames(data), Cond=noRep(colnames(data)))
    new("AnnotatedDataFrame", data=pData)
}

runTest <- function(data, test)
{
    case <- data[, noRep(colnames(data)) == test$Mortal[[1]]]   # Mortal
    cont <- data[, noRep(colnames(data)) == test$Immortal[[1]]] # Immortal
    
    if (ncol(case) != 2 && ncol(case) != 3) { return(NULL) }
    if (ncol(cont) != 2 && ncol(cont) != 3) { return(NULL) }    
    
    stopifnot(length(unique(noRep(colnames(case)))) == 1)
    stopifnot(length(unique(noRep(colnames(cont)))) == 1)    
    
    # Always keep controls before cases
    newData <- cbind(cont, case)
    
    #
    # http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html
    #
    
    if (ncol(cont) == 2) { design <- model.matrix(~factor(c(1,1,2,2,2))) }
    if (ncol(cont) == 3) { design <- model.matrix(~factor(c(1,1,1,2,2,2))) }
    
    colnames(design) <- c("Intercept", "Diff")
    eb.fit(newData, design)
}

checkMissing <- function()
{
    checkTest <- function(data, test)
    {
        quantile(data[!is.na(rowSums(data)),])
        
        print("---------------------------------")
        print(test$Name)
        case <- data[, noRep(colnames(data)) == test$Mortal[[1]]]   # Mortal
        cont <- data[, noRep(colnames(data)) == test$Immortal[[1]]] # Immortal

        if (ncol(cont) == 0) { print("Missing all replicates in immortal") }
        if (ncol(case) == 0) { print("Missing all replicates in mortal")   }
        
        mCont <- rowSums(is.na(cont))
        mCase <- rowSums(is.na(case))        
        
        print(sum(mCont == 1))
        print(sum(mCont == 2))
        print(sum(mCont == 3))
        print(sum(mCase == 1))
        print(sum(mCase == 2))        
        print(sum(mCase == 3))
    }

    tests <- getTests()
    tests$Name <- paste(tests$Mortal, tests$Immortal, sep=" - ")
    data  <- read.table("SWATH/data2.tsv", row.names=1, header=TRUE, sep='\t')
    print(dim(data)) # 5854 x 110
    for (i in 1:nrow(tests)) { checkTest(data, tests[i,]) }
}

#min(data, na.rm=TRUE)
#quantile(data, na.rm=TRUE)
#print(data[is.na(rowSums(data)), is.na(colSums(data))])
#checkMissing()
