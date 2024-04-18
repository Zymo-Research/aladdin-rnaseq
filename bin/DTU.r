#!/usr/bin/env Rscript

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    usage_str <- paste("Usage: DTU.r <Salmon results directory> <design CSV> <FDR cutoff>",
    "<min. proportion of samples where transcript has >=10 reads> <min. proportion of samples where transcript accounts for >=10% expression of the gene>", 
    "<min. proportion of samples where gene has >=10 reads> <comparisons CSV(optional)")
    stop(usage_str, call.=FALSE)
}
input <- args[1]
design_csv <- args[2]
fdr_cutoff <- as.numeric(args[3])
prop_filter_transcript_counts <- as.numeric(args[4])
prop_filter_transcript_props <- as.numeric(args[5])
prop_filter_gene_counts <- as.numeric(args[6])

# Load / install packages
if (!require("DEXSeq")) {
    install.packages("BiocManager")
    BiocManager::install("DEXSeq")
    library("DEXSeq")
}
if (!require("tximport")) {
    install.packages("BiocManager")
    BiocManager::install("tximport")
    library("tximport")
}
if (!require("DRIMSeq")) {
    install.packages("BiocManager")
    BiocManager::install("DRIMSeq")
    library("DRIMSeq")
}
if (!require("stageR")) {
    install.packages("BiocManager")
    BiocManager::install("stageR")
    library("stageR")
}
if (!require("openxlsx")) {
    install.packages("openxlsx")
    library("openxlsx")
}

# Read the sample grouping
grouping <- read.csv(design_csv, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, na.strings='')
# Read the tx2gene file
tx2gene <- read.csv('tx2gene.tsv', header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE, na.strings='')
# Extract gene names if available, if not, use gene_id
if ('gene_name' %in% colnames(tx2gene)) {
    gene_name <- tx2gene[,c('gene_id','gene_name')]
    gene_name <- gene_name[!duplicated(gene_name),]
    gene_name <- data.frame(gene_name$gene_name, row.names=gene_name$gene_id, check.names=FALSE, stringsAsFactors=FALSE)
} else {
    gene_name <- unique(tx2gene$gene_id)
    gene_name <- data.frame(gene_name, row.names=gene_name, check.names=FALSE, stringsAsFactors=FALSE)
}
# List samples in Salmon results dir
samples <- list.dirs(input, recursive=FALSE, full.names=FALSE)
# Get a list of samples that are shared by the input and design file
shared_samples <- intersect(grouping$sample, samples)
# Import Salmon results
files <- file.path(input, shared_samples, 'quant.sf')
names(files) <- shared_samples
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut=TRUE, countsFromAbundance = "dtuScaledTPM")

# Clean up the grouping file so it only contains samples that are shared by the input and design file
grouping <- grouping[grouping$sample %in% shared_samples,]
# If group label absent, use sample label
grouping$group <- ifelse(is.na(grouping$group), grouping$sample, grouping$group)
# Filter groups so that all groups have at least 2 replicates
no.replicates <- table(grouping$group)
groups <- names(no.replicates)[no.replicates > 1]
# Less than 2 groups with replicates, exit gracefully.
if (length(groups) < 2) {
    message('There are fewer than two groups with replicates. Cannot carry out DTU analysis.')
    quit(save='no', status=0)
}
# Read requested group comparisons, or run every pairwise group comparisons
if (length(args) > 6) {
    comparisons_file <- args[7]
    comparisons <- read.csv(comparisons_file, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, na.strings='')
    comparisons <- t(comparisons)
    rownames(comparisons) <- c()
} else {
    comparisons <- combn(groups, 2)
}

# Run DTU on every pair of groups
for (i in 1:ncol(comparisons)) {
    g1 <- comparisons[1,i]
    g2 <- comparisons[2,i]
    if (!(g1 %in% groups) || !(g2 %in% groups)) {
        message(paste(g1, "or", g2, "did not have replicates. Skipping comparison..."))
        next
    }
    # Subset data based on the two groups
    grouping_subset <- grouping[grouping$group %in% c(g1,g2),]
    sample_subset <- grouping[grouping$group %in% c(g1,g2),'sample']
    cts <- txi$counts[,sample_subset]
    # Filter out transcripts with zero counts
    cts <- cts[rowSums(cts)>0, ]
    tx2gene_subset <- tx2gene[match(row.names(cts), tx2gene$transcript_id),]
    # Construct DRIMSeq object
    counts <- data.frame(gene_id=tx2gene_subset$gene_id, feature_id=tx2gene_subset$transcript_id, cts)
    samps <- data.frame(sample_id=grouping_subset$sample, group=grouping_subset$group)
    d <- dmDSdata(counts=counts, samples=samps)
    # Filter out transcripts with low read counts or low proportion of counts within the same gene
    no.samples <- length(sample_subset)
    no.samples.filter_transcript_counts <- ceiling(no.samples*prop_filter_transcript_counts)
    no.samples.filter_transcript_props <- ceiling(no.samples*prop_filter_transcript_props)
    no.samples.filter_gene_counts <- ceiling(no.samples*prop_filter_gene_counts)
    d <- dmFilter(d, 
        min_samps_feature_expr=no.samples.filter_transcript_counts, min_feature_expr=10, 
        min_samps_feature_prop=no.samples.filter_transcript_props, min_feature_prop=0.1, 
        min_samps_gene_expr=no.samples.filter_gene_counts, min_gene_expr=10)

    # Construct DEXSeq object
    sample.data <- DRIMSeq::samples(d)
    count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
    dxd <- DEXSeqDataSet(countData=count.data,
                        sampleData=sample.data,
                        design=~sample + exon + group:exon,
                        featureID=counts(d)$feature_id,
                        groupID=counts(d)$gene_id)

    # Run DEXseq analysis
    dxd <- estimateSizeFactors(dxd)
    dxd <- estimateDispersions(dxd, quiet=TRUE)
    dxd <- testForDEU(dxd, reducedModel=~sample + exon)
    dxd <- estimateExonFoldChanges(dxd, fitExpToVar="group", denominator=g2)
    dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
    qval <- perGeneQValue(dxr)

    # Run stageR for better OFDR control
    pConfirmation <- matrix(dxr$pvalue, ncol=1)
    dimnames(pConfirmation) <- list(dxr$featureID, "transcript")
    stageRObj <- stageRTx(pScreen=qval, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene_subset)
    stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=fdr_cutoff)
    suppressWarnings({ dex.padj <- getAdjustedPValues(stageRObj, order=TRUE, onlySignificantGenes=FALSE) })

    # Combine results and output
    colnames(dex.padj)[c(3,4)] <- c('screening_padj','confirmation_padj')
    dxr <- as.data.frame(dxr[, !names(dxr) %in% c('genomicData','countData')])
    results <- merge(dxr, dex.padj, by.x=c('groupID','featureID'), by.y=c('geneID','txID'), all=TRUE)
    results <- merge(gene_name, results, by.x=0, by.y=1, all.y=TRUE)
    colnames(results)[1:3] <- c('Gene ID', 'Gene Name', 'Transcript ID')
    results <- results[order(results$screening_padj, results$'Gene ID', results$confirmation_padj),]
    write.table(results, paste0(g1,'_vs_',g2,'_DTU_analysis_DEXSeq_results.tsv'), sep="\t", quote=FALSE, row.names=FALSE)
    write.xlsx(results, paste0(g1,'_vs_',g2,'_DTU_analysis_DEXSeq_results.xlsx'))
}