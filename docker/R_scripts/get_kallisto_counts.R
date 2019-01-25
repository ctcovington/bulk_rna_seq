args = commandArgs(trailingOnly=TRUE)
transcript_counts_dir <- args[1]
organism <- args[2]
samples_described_path <- args[3]
samples_compared_path <- args[4]
output_dir <- args[5]
stopifnot(organism %in% c('human', 'mouse'))

library(rhdf5)
library(tximport)
library(biomaRt)
library(readr)
library(data.table)
library(sleuth)
library(stringr)
library(jsonlite)

if (organism == 'human') {
    # download transcript to gene database
    library(EnsDb.Hsapiens.v86)
    edb <- EnsDb.Hsapiens.v86

    # download ensembl database
    ensembl <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
} else if (organism == 'mouse') {
    # download transcript to gene database
    library(EnsDb.Mmusculus.v79)
    edb <- EnsDb.Mmusculus.v79

    # download ensembl database
    ensembl <- useDataset('mmusculus_gene_ensembl', useMart('ensembl'))
}

# ####################################################################
# ### sleuth version                                               ###
# ####################################################################
# # run Sleuth if and only if samples_described and samples_compared are not empty
# if (file.info(samples_described_path)$size != 0 & file.info(samples_compared_path)$size != 0) {
#     # create transcript to symbol mapping
#     if (organism == 'human') {
#         mart <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
#                                  dataset = 'hsapiens_gene_ensembl',
#                                  host = 'ensembl.org')
#         t2g <- getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol'), mart = mart)
#         t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, gene = hgnc_symbol)
#
#     } else if (organism == 'mouse') {
#         mart <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
#                                  dataset = 'mmusculus_gene_ensembl',
#                                  host = 'ensembl.org')
#         t2g <- getBM(attributes = c('ensembl_transcript_id', 'mgi_symbol'), mart = mart)
#         t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, gene = mgi_symbol)
#     }
#     # NOTE: some transcripts map to multiple genes -- so here I combine the gene names for those transcripts
#     t2g <- aggregate(gene ~ target_id, data = t2g, paste, sep = '_')
#
#     # get list of sample IDs
#     sample_id <- dir(transcript_counts_dir)
#
#     # create list of paths to each directory
#     kal_dirs <- file.path(transcript_counts_dir, sample_id)
#
#     # load sample information
#     s2c <- data.table(read.table(samples_described_path))
#     s2c[, path := kal_dirs]
#     setnames(s2c, names(s2c), c('condition', 'sample', 'path'))
#
#     # iterate over samples and remove those with < 30% alignment
#     for (directory in s2c[, path]) {
#         # extract sample name
#         sample_name <- str_replace(basename(directory), '_S[0-9]_001', '')
#
#         # load run info
#         qc <- fromJSON(file.path(directory, 'run_info.json'))
#
#         # get percentage pseudoaligned
#         percentage_aligned <- qc$p_pseudoaligned
#
#         # remove samples with < 30% alignment
#         if (percentage_aligned < 30) {
#             s2c <- s2c[sample != sample_name]
#         }
#     }
#
#     # generate list of differential expressioncomparisons
#     comparisons <- fread(samples_compared_path, header = FALSE, sep = '\t')
#
#     # for each comparison; initialize sleuth object, set up models, run differential expression, and output table
#     for (comparison in comparisons$V1) {
#         cat(sprintf('Performing Sleuth differential expression for %s\n\n', comparison))
#         # identify relevant comparison and samples
#         conditions <- str_split(comparison, '_vs_')[[1]]
#         samples_a <- which(s2c$condition == conditions[1])
#         samples_b <- which(s2c$condition == conditions[2])
#
#         # initialize sleuth object by identifying samples of interest and aggregating counts at the gene level
#         so <- sleuth_prep(s2c[c(samples_a, samples_b)], target_mapping = t2g, aggregation_column = 'gene', 'gene_mode' = TRUE, extra_bootstrap_summary = TRUE)
#
#         # fit models
#         so <- sleuth_fit(so, ~ condition, 'full') # fit full model
#         so <- sleuth_fit(so, ~ 1, 'reduced') # fit reduced model
#         so <- sleuth_lrt(so, 'reduced', 'full') # test difference in models
#
#         # create differential expression table
#         sleuth_table <- data.table(sleuth_results(so, 'reduced:full', 'lrt', pval_aggregate = FALSE, show_all = FALSE))
#
#         # clean up differential expression table
#         setnames(sleuth_table, 'target_id', 'gene')
#         sleuth_table <- sleuth_table[gene != ''] # remove unnamed genes
#
#         # print sleuth_table to file
#         fwrite(sleuth_table, file.path(output_dir, sprintf('%s_diff_exp.tsv', comparison)), sep = '\t')
#
#         cat('\n')
#     }
# }


####################################################################
### edgeR version                                                ###
####################################################################
# get list of sample IDs
sample_id <- dir(transcript_counts_dir)

# create list of paths to each directory
kal_dirs <- file.path(transcript_counts_dir, sample_id)

#######################################
### create transcript to gene table ###
#######################################
# create table mapping between transcripts and ensembl gene ID
tx2gene <- data.table(transcripts(edb, return.type = 'data.frame'))
tx2gene <- tx2gene[, c('tx_id', 'gene_id')]
setnames(tx2gene, names(tx2gene), c('TXNAME', 'GENEID'))

#############################
## create counts for edgeR ##
#############################
abundance_files <- file.path(kal_dirs, 'abundance.h5')
names(abundance_files) <- sample_id
kallisto_object <- tximport(abundance_files, type = 'kallisto', tx2gene = tx2gene, ignoreTxVersion = TRUE)
gene_names <- row.names(kallisto_object$counts)
kallisto_counts <- data.table(kallisto_object$counts)
kallisto_counts[, gene := gene_names]

# replace ensembl gene ID with HGNC/MGI symbol
# NOTE: hgnc/mgi symbol will still be under the name 'GENEID' because that is what tximport expects
if (organism == 'human') {
    mart <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
                             dataset = 'hsapiens_gene_ensembl',
                             host = 'ensembl.org')
    gene_map <- data.table(getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart = mart))
} else if (organism == 'mouse') {
    mart <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
                             dataset = 'mmusculus_gene_ensembl',
                             host = 'ensembl.org')
    gene_map <- data.table(getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'), mart = mart))
}
setnames(gene_map, 'ensembl_gene_id', 'gene')
kallisto_counts <- merge(x = kallisto_counts, y = gene_map, all.x = TRUE, by = c('gene'))
if (organism == 'human') {
    kallisto_counts[, gene := hgnc_symbol]
    kallisto_counts[, hgnc_symbol := NULL]
} else if (organism == 'mouse') {
    kallisto_counts[, gene := mgi_symbol]
    kallisto_counts[, mgi_symbol := NULL]
}

# sum over rows with same hgnc/mgi symbol
count_cols <- setdiff(names(kallisto_counts), 'gene')
kallisto_counts[, (count_cols) := lapply(.SD, sum), by = gene, .SDcols = count_cols]
kallisto_counts <- unique(kallisto_counts)

# remove missing genes
kallisto_counts <- kallisto_counts[gene != '' & !is.na(gene)]

# sort alphabetically by gene
kallisto_counts <- kallisto_counts[order(gene)]

# round counts and clean count names
count_names <- setdiff(names(kallisto_counts), 'gene')
kallisto_counts[, (count_names) := lapply(.SD, round), .SDcols = count_names]
setnames(kallisto_counts, count_names, str_replace(count_names, '_001', ''))

# write to file
fwrite(kallisto_counts, file.path(output_dir, 'kallisto_counts.csv'))
