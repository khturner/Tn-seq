Args <- commandArgs(TRUE)
#Args <- c("Control", "ControlReps", "Test", "TestReps", "AssemblyPrefix", "OutputPrefix")
sites <- read.table(paste(Args[6], "-sites.txt", sep=""), header=TRUE)
sites.ordered <- sites[order(sites$Position),]

# LOESS smooth data
for (i in 2:(length(sites.ordered))) { 
	counts.loess <- loess(sites.ordered[[i]] ~ sites.ordered$Position, span=1, data.frame(x=sites.ordered$Position, y=sites.ordered[[i]]), control=loess.control(statistics=c("approximate"),trace.hat=c("approximate")))
	counts.predict <- predict(counts.loess, data.frame(x=sites.ordered$Position))
	counts.ratio <- counts.predict/median(counts.predict)
	sites.ordered[[i]] <- sites.ordered[[i]]/counts.ratio
}
countsmatrix <- matrix(sites.ordered[[2]])
for (i in 3:(length(sites.ordered))) {
	newmatrix <- cbind(countsmatrix, sites.ordered[[i]])
	countsmatrix <- newmatrix
}
rm(newmatrix, i)

gff <- read.delim(file=paste(Args[5], ".trunc.gff", sep=""), sep="\t")
colnames(gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "att", "x")
gff <- tail(gff, n=-5)
library(DESeq)
sitescds <- newCountDataSet(round(countsmatrix), c(rep(Args[1], as.numeric(Args[2])), rep(Args[3], as.numeric(Args[4]))))
sitescds <- estimateSizeFactors(sitescds)

counts.norm <- cbind(sites.ordered[,1],t(t(counts(sitescds))/sizeFactors(sitescds)))
controlreps <- 0
testreps <- 0
numsites <- matrix(0, nrow=length(gff[,1]), ncol=length(counts.norm[1,])-1)
for (c in 2:length(counts.norm[1,])) {
	for (i in 1:length(gff[,1])) {
		gff[i,c+8] <- 1
	}
	if (controlreps < Args[2]) {
		controlreps <- controlreps + 1
		colnames(gff)[c+8] <- paste(Args[1], controlreps)
	}
	else {
		testreps <- testreps + 1
		colnames(gff)[c+8] <- paste(Args[3], testreps)
	}
}

line <- 1
for (r in 1:length(counts.norm[,1])) {
	pos <- counts.norm[r,1]
	while (line <= length(gff$start) & gff$start[line] <= pos) {
		if (gff$end[line] >= pos) {
			for (c in 2:length(counts.norm[1,])) {
				gff[line,c+8] <- gff[line,c+8] + counts.norm[r,c]
				if (counts.norm[r,c] > 0) {
					numsites[line,c-1] <- numsites[line,c-1] + 1
				}
			}
		}
		line <- line + 1
	}
	while (line > 1 & gff$start[line] > pos) {
		line <- line - 1
	}
}

genecounts <- gff[,10]
for (i in 11:(length(gff[1,]))) {
	newmatrix <- cbind(genecounts, gff[,i])
	genecounts <- newmatrix
}
genes <- read.delim(file=paste(Args[5], ".gene.products.kegg.txt", sep=""), sep="\t", header=TRUE)
rownames(genecounts) <- genes$Locus
colnames(genecounts) <- colnames(gff)[10:length(gff)]
colnames(counts.norm) <- c("Position", colnames(gff)[10:length(gff)])
colnames(numsites) <- colnames(gff)[10:length(gff)]
genescds <- newCountDataSet(round(genecounts), c(rep(Args[1], as.numeric(Args[2])), rep(Args[3], as.numeric(Args[4]))))
genescds$sizeFactor <- rep(1, length(genecounts[1,]))
#genescds <- estimateDispersions(genescds, fitType="local", sharingMode="fit-only") #use this if estimateDispersions fails
genescds <- estimateDispersions(genescds)
res <- nbinomTest(genescds, Args[1], Args[3])
colnames(res)[3] <- paste(Args[1], "Mean", sep="")
colnames(res)[4] <- paste(Args[3], "Mean", sep="")
res <- cbind(res, genes[,2:5], numsites)

write.table(counts.norm, file=paste(Args[6], ".sites.counts.tsv", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
write.table(res, file=paste(Args[6], ".DESeq.tsv", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
