# Defined as a function, return value is table of differential fitness results with essentiality calls (source into R/RStudio for interactive analyses)

TnSeqDESeqEssential <- function(ctrl_pfx, ctrl_reps, gff_pfx, out_pfx, to_trim, num_expected, in_files) {
	# Read in sites files
	library(dplyr)
	sites <- data.frame(Pos=c(0)) %>% tbl_df
	for (i in 1:length(in_files)) {
		newsites <- read.table(paste(paste(in_files[i], in_files[i], sep="/"), "sites.txt", sep="-")) %>% tbl_df
		colnames(newsites) <- c(paste("V", i, sep=""), "Pos")
		newsites <- tail(newsites, n=-to_trim) %>% arrange(Pos)
		sites <- merge(sites, newsites, all=T) %>% tbl_df
	}
	sites <- tail(sites, n=-1)
	sites[is.na(sites)] <- 0

	# OPTIONAL - perform site filtering. Example: only consider sites identified in both of 2 replicates
	#sites <- sites %>% mutate(numreps = (V1 > 0) + (V2 > 0)) %>% filter(numreps == 2)
	#sites <- sites[-4]
	
	# LOESS smooth data
	for (i in 2:(length(sites))) { 
		counts.loess <- loess(sites[,i] ~ sites$Pos, span=1, data.frame(x=sites$Pos, y=sites[,i]), control=loess.control(statistics=c("approximate"),trace.hat=c("approximate")))
		counts.predict <- predict(counts.loess, data.frame(x=sites$Pos))
		counts.ratio <- counts.predict/median(counts.predict)
	    sites[,i] <- sites[,i]/counts.ratio
	}

	# Normalize data by reads/site
	library(DESeq)
	sitescds <- sites[,2:length(sites)] %>% round %>% newCountDataSet(c(rep(ctrl_pfx, ctrl_reps)))
	sitescds <- estimateSizeFactors(sitescds)
	counts.norm <- counts(sitescds, normalized=T)
	rownames(counts.norm) <- sites$Pos
	
	# Initialize the list of genes, determine genome length
	gff <- read.delim(file=paste(gff_pfx, ".trunc.gff", sep=""), sep="\t", fill=TRUE) %>% tbl_df
	colnames(gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "att", "x")
	genomelength <- as.numeric(strsplit(as.character(gff[3,1]), " ")[[1]][4])
	gff <- tail(gff, n=-5)
	gff <- gff[(gff$feature=="gene"),]

	# Generate pseudo-datasets with the same number of insertion sites and total reads mapping to those sites, randomly distributed across the genome
	print("Generating pseudo-datasets")
	counts.df <- data.frame(counts.norm)
	counts.df$Pos <- as.numeric(rownames(counts.df))
	numreads <- sum(counts.norm)/ctrl_reps
	numsites <- length(which(counts.norm>0))/ctrl_reps
	for (i in 1:num_expected) {
		expected <- data.frame(Pos=sample(1:genomelength, numsites), Exp=sample(sites$V1, numsites)) %>% tbl_df %>% arrange(Pos)
		colnames(expected)[2] <- paste("Expected", i, sep=".")
		counts.df <- merge(counts.df, expected, by="Pos", all=T) %>% tbl_df %>% arrange(Pos)
		counts.df[is.na(counts.df)] <- 0
	}
	rownames(counts.df) <- counts.df$Pos
	counts.norm <- as.matrix(counts.df[,(2:length(counts.df))])

	# Initialize the lists of read counts per gene and number of independent Tn sites per gene
	controlreps <- 0
	expreps <- 0
	for (c in 1:length(counts.norm[1,])) {
		gff[,c+9] <- rep(1,length(gff[,1]))
		if (controlreps < ctrl_reps) {
			controlreps <- controlreps + 1
			colnames(gff)[c+9] <- paste(ctrl_pfx, controlreps, sep=".")
		}
		else {
			expreps <- expreps + 1
			colnames(gff)[c+9] <- paste("Expected", expreps, sep=".")
		}
	}

	# Output gene boundaries and read counts per Tn site for Perl binning script
	print("Binning read counts by gene boundaries")
	boundariesfile <- paste(out_pfx, ".boundaries.tsv", sep="")
	sitecountsfile <- paste(out_pfx, ".sitecounts.tsv", sep="")
	write.table(gff[,c(4,5, 10:length(gff))], boundariesfile, quote=FALSE, sep="\t", row.names=F)
	write.table(counts.df, sitecountsfile, quote=FALSE, sep="\t", row.names=F)
	system(paste("TnGeneBin.pl", boundariesfile, sitecountsfile))
	genecounts <- read.table(paste(boundariesfile, "out", sep="."), header=T)[,-c(1,2)]
	numsites <- read.table(paste(boundariesfile, "numsites.out", sep="."), header=T)[,-c(1,2)]
	system(paste("rm", boundariesfile,
		paste(boundariesfile, "out", sep="."),
		paste(boundariesfile, "numsites.out", sep=".")))

	# Uncomment this section if you have a kegg annotation description file of the genes and their products
	genes <- read.delim(file=paste(gff_pfx, ".gene.products.kegg.txt", sep=""), sep="\t", header=TRUE)
	rownames(genecounts) <- genes$Locus
	
	
	# Uncomment this section if you DO NOT have a kegg annotation description file of the genes and their products
	#genes <- matrix("", length(gff[,1]), 2)
	#for (i in 1:length(gff[,1])) {
	#	genes[i,1] <- strsplit(grep("locus_tag",strsplit(as.character(gff$att[i]),";")[[1]], value=T),"=")[[1]][2]
	#	genes[i,2] <- strsplit(grep("product",strsplit(as.character(gff$att[i]),";")[[1]], value=T),"=")[[1]][2]
	#}
	write.table(genecounts, paste(out_pfx, ".genecounts.tsv", sep=""), quote=FALSE, sep="\t")

	# Perform differential fitness analysis
	colnames(numsites) <- colnames(gff)[10:length(gff)]
	numsitesout <- data.frame(numsites[,(1:ctrl_reps)])
	numsitesout[,ctrl_reps+1] <- rowMeans(numsites[,-(1:ctrl_reps)])
	colnames(numsitesout)[ctrl_reps+1] <- "Expected"
	genescds <- newCountDataSet(round(genecounts), c(rep(ctrl_pfx, ctrl_reps), rep("Expected", num_expected)))
	genescds$sizeFactor <- rep(1, length(genecounts[1,])) # This is manually set as 1 because we normalized by site above
	genescds <- estimateDispersions(genescds)
	res <- nbinomTest(genescds, "Expected", ctrl_pfx) %>% tbl_df
	colnames(res)[4] <- paste(ctrl_pfx, "Mean", sep="")
	colnames(res)[3] <- "ExpectedMean"
	out <- cbind(res, genes[,2:5], numsitesout) %>% tbl_df # Uncomment if you have a kegg annotation
	#out <- cbind(res, genes[,2:3], numsitesout) %>% tbl_df # Uncomment if you do not have a kegg annotation

	# Perform bimodal clustering and essentiality calling and output results
	library(mclust)
	fit <- Mclust(out$log2FoldChange, G=1:2)
	category <- rep("",length(out$id))
	for (i in 1:length(out$id)) {
		if (fit$classification[i] == 1 & out$log2FoldChange[i] < 0) {
			category[i] <- "Reduced"
		}
		else {
			category[i] <- "Unchanged"
		}
	}
	fit$uncertainty[which(out$log2FoldChange > 0)] <- 0
	essentiality <- cbind(category, fit$uncertainty)
	colnames(essentiality) <- c("Essentiality", "Uncertainty")
	out <- cbind(out, essentiality) %>% tbl_df
	write.table(out, file=paste(out_pfx, ".DESeq.tsv", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
	return(out)
}
Args <- commandArgs(TRUE)
TnSeqDESeqEssential(Args[1], as.numeric(Args[2]), Args[3], Args[4], as.numeric(Args[5]), as.numeric(Args[6]), Args[-(1:6)])
