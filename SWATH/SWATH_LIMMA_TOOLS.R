library(limma)
library(ggplot2)

eb.fit <- function(dat, design, mortal, immortal, saveTo)
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
    
    data <- data.frame(name=row.names(dat), lqm=-log2(q.mod), qm=q.mod, logFC=logFC)
    data$color <- "black"

    if (nrow(data[data$qm <= 0.05,]) > 0) { data[data$qm <= 0.05,]$color <- "red"    }
    if (nrow(data[abs(logFC) >= 3,]) > 0) { data[abs(logFC) >= 3,]$color <- "orange" }
    if (nrow(data[data$qm <= 0.05 & abs(logFC) >= 3,]) > 0) { data[data$qm <= 0.05 & abs(logFC) >= 3,]$color <- "green" }
    data$color <- as.factor(data$color)

    png(paste("~/Desktop/", saveTo, sep=""))
    p <- ggplot(data, aes(y=lqm, x=logFC))
    p <- p + geom_point(aes(col=color))
    p <- p + scale_colour_manual(values = c(levels(data$color)))
    p <- p + geom_vline(xintercept=4, color="black", size=0.2)
    p <- p + geom_vline(xintercept=-4, color="black", size=0.2)
    p <- p + theme_bw()
    p <- p + xlab("Log fold changes")
    p <- p + ylab("-log2 moderated q-value")
    p <- p + ggtitle(paste("Mortal: ", mortal, "    Immortal: ", immortal, sep=""))
    p <- p + theme(legend.position="none")
    print(p)
    dev.off()

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

runTest <- function(data, test, saveTo)
{
    norm <- data[, noRep(colnames(data)) == test$Mortal[[1]]]   # Mortal (normal controls)
    tumo <- data[, noRep(colnames(data)) == test$Immortal[[1]]] # Immortal (tumor sample)
    
    stopifnot(sum(is.na(norm)) == 0)
    stopifnot(sum(is.na(tumo)) == 0)    

    if (ncol(norm) == 0) { return(NULL) }
    if (ncol(tumo) == 0) { return(NULL) }    
    
    stopifnot(ncol(norm) == 2 || ncol(norm) == 3)
    stopifnot(ncol(tumo) == 2 || ncol(tumo) == 3)    

    stopifnot(length(unique(noRep(colnames(norm)))) == 1)
    stopifnot(length(unique(noRep(colnames(tumo)))) == 1)
    
    design <- c(rep(-1, ncol(norm)), rep(1, ncol(tumo)))
    design <- model.matrix(~factor(design))
    colnames(design) <- c("Intercept", "Diff")

    r <- eb.fit(cbind(norm, tumo), design, test$Mortal[[1]], test$Immortal[[1]], saveTo)
    
    acc  <- function(x) { sapply(strsplit(x, split="|", fixed=TRUE) , "[[", 2) }
    gene <- function(x) { sapply(strsplit(x, split="|", fixed=TRUE) , "[[", 3) }    
    
    d <- data.frame(row.names = row.names(r),
                    access    = acc(row.names(r)),
                    gene      = gene(row.names(r)),
                    mortalM   = rowMeans(norm),
                    immortM   = rowMeans(tumo),
                    logFC     = r$logFC,
                    tstats    = r$t.ord,
                    tmod      = r$t.mod,
                    pval      = r$p.ord,
                    pmod      = r$p.mod,
                    qval      = r$q.ord,
                    q.mod     = r$q.mod)

    d$mortal1   = norm[,1]
    d$mortal2   = norm[,2]
    if (ncol(norm) == 2) { d$mortal3   <- rep(NaN, nrow(data)) }
    if (ncol(norm) == 3) { d$mortal3   <- norm[,3] }

    d$immortal1 = tumo[,1]
    d$immortal2 = tumo[,2]
    if (ncol(tumo) == 2) { d$immortal3 <- rep(NaN, nrow(data)) }
    if (ncol(tumo) == 3) { d$immortal3 <- tumo[,3] }
    
    d$MortalName   <- test$Mortal[[1]]
    d$ImmortalName <- test$Immortal[[1]]
    d
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