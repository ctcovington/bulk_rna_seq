library(limma)
library(edgeR)
library(stringr)
library(data.table)

filter_cpm_cutoff <- function(counts, cpm.cutoff=100, min.samples=1) {
  # keep genes that pass cpm.cutoff in at least min.samples
  keep.rows = (rowSums(cpm(counts) >= cpm.cutoff) >= min.samples)
  return (counts[keep.rows,])
}

filter_cpm_rank <- function(counts, keep.genes=8000, method='var') {
  # keep genes that rank highest in some measure of cross-sample CPM
  rank.metrics = c('var','max','sum')
  keep.genes = min(keep.genes, nrow(counts))
  if (method %in% rank.metrics) {
    ranks = apply(cpm(counts), 1, method)
    ranks = sort(ranks, decreasing=TRUE)
    genes = names(ranks)[1:keep.genes]
    keep.rows = (rownames(counts) %in% genes)
    return (counts[keep.rows,])
  } else {
    print(paste('Unrecognized method', method, '- skipping filter'))
    return (counts)
  }
}

filter_zeros <- function(counts) {
    return (counts[which(rowSums(counts) > 0),])
}

get_geomean_libsizes <- function(counts, targets=NULL) {
    if (is.null(targets)) {
    targets = c('Actb', 'B2m', 'Gapdh', 'Tbca', 'Ubc', 'Ywhaz', 'Zyg11b')
    }
    targets = targets[which(targets %in% rownames(counts))]
    geom.libsizes = geom_mean_cols(counts[targets,]/colSums(counts))
    # divide by max, so that the largest library size is 1
    # thus when dividing by these, each sample's new 'pseudocounts'
    # should be in the same order of magnitude as original counts
    geom.libsizes = geom.libsizes / max(geom.libsizes)

    # TODO this only gives per-sample, not per-group norms
    return (geom.libsizes)
}

geom_mean <- function(x) {
  # only works for counts (non-negative)
  # likely to overflow for large inputs
  x.eff = x[which(x > 0)]
  p = prod(x.eff)
  n = length(x.eff)
  return (p^(1/n))
}

pseudo_geom_mean <- function(x) {
  return (exp(mean(log(x+1))))
}

geom_mean_cols <- function(matrix) {
  gmc = rep(1,times=ncol(matrix))
  for (ii in 1:ncol(matrix)) {
    gmc[ii] = geom_mean(matrix[,ii])
  }
  return (gmc)
}

filter_design <- function(design, subset) {
  # filter design matrix, keeping only relevant DOF
  design.new = design[subset,]
  keep.dof = !apply(design.new, 2, function(x){all(diff(x)==0)})
  keep.dof[1] = TRUE # always keep intercept
  return (design.new[,keep.dof])
}

run_edger_glm <- function(counts, design, keep.samples=NULL, calc.norm=TRUE, out.file=NULL) {
  if (!is.null(keep.samples)) {
    design.new = filter_design(design, keep.samples)
    dge = DGEList(counts[,keep.samples])
  } else {
    design.new = design
    dge = DGEList(counts)
  }
  if (calc.norm) { dge = calcNormFactors(dge) }

  # check to see if and groups have only one replicate -- if so, perform modified modeling
  if (any(colSums(design)==1)) {
      dge = estimateGLMCommonDisp(dge, method = 'deviance', robust = TRUE, subset = NULL)
      dge = estimateGLMTrendedDisp(dge, method = 'auto')
      dge = estimateGLMTagwiseDisp(dge)
  } else {
      dge = estimateGLMCommonDisp(dge, design.new)
      dge = estimateGLMTrendedDisp(dge, design.new)
      dge = estimateGLMTagwiseDisp(dge, design.new)
  }
  fit = glmFit(dge, design.new)
  lrt = glmLRT(fit)
  if (!is.null(out.file)) {
    ranked = topTags(lrt, n=NULL)
    write.table(ranked, file=out.file, sep='\t', quote=F, row.names=T)
  }
  return (lrt)
}

run_edger_classic <- function(counts, a.cols, b.cols, calc.norm=TRUE, out.file=NULL) {
  grouping = factor(c(rep(0,length(a.cols)), rep(1,length(b.cols))))
  dge = DGEList(counts[,c(a.cols,b.cols)], group=grouping)
  if (calc.norm) { dge = calcNormFactors(dge) }
  dge = estimateCommonDisp(dge)
  dge = estimateTagwiseDisp(dge)
  et = exactTest(dge)
  if (!is.null(out.file)) {
    ranked = topTags(et, n=NULL)
    write.table(ranked, file=out.file, sep='\t', quote=F, row.names=T)
  }
  return (et)
}

run_single_glm_analysis <- function(counts, design.file, primary, secondaries, out.file) {
  # NB: secondaries can be empty/NULL, but still a required arg
  design.table = read.delim(design.file, row.names=1)
  design.formula = as.formula(paste('~',paste(c(secondaries, primary), collapse='+'), sep=''))
  design.model = model.matrix(design.formula, data=design.table)
  design.samples = rownames(design.table)
  counts.subset = counts[,which(colnames(counts) %in% design.samples)]
  results = run_edger_glm(counts.subset, design.model, out.file=out.file)
}

run_full_glm_analysis <- function(counts.file, design.listing.file, design.dir, out.dir, keep.genes=8000) {
  counts = read.delim(counts.file, row.names=1)
  counts = filter_cpm_rank(counts, keep.genes=keep.genes, method='var')
  design.listing = read.delim(design.listing.file, header=FALSE)

  for (ii in 1:nrow(design.listing)) {
    basename = as.character(design.listing[ii,1])
    primary = as.character(design.listing[ii,2])
    if (length(design.listing[ii,]) >= 3) {
      secondaries = unlist(strsplit(as.character(design.listing[ii,3]),','))
    } else {
      secondaries = NULL
    }
    design.file = file.path(design.dir, basename)
    out.file = file.path(out.dir, paste(basename,'edgeR','tsv',sep='.'))
    run_single_analysis(counts, design.file, primary, secondaries, out.file)
  }

}

get_simple_design <- function(samples.described.file) {
  samples = read.table(samples.described.file, row.names=2, header=FALSE, check.names=F)
  design = model.matrix(~0+samples$V1)
  colnames(design) = make.names(levels(samples$V1))
  rownames(design) = rownames(samples)
  #  print(design)
  return (design)

}

get_complex_design <- function(samples.described.file) {
  samples = read.delim(samples.described.file, row.names=1, header=FALSE, check.names=F)
  design = model.matrix(~0+samples$V1+samples$V2)
  colnames(design) = make.names(levels(samples$V1))
  rownames(design) = rownames(samples)
  #  print(design)
  return (design)

}

get_simple_contrasts <- function(samples.compared.file, design) {
  comparison.names = readLines(samples.compared.file)
  comparison.names.sanitized = make.names(comparison.names)
  comparison.pairs = gsub('_vs_', ' ', comparison.names.sanitized)
  comparison.pairs = read.table(textConnection(comparison.pairs))
  # our current convention is to map A_vs_B to the model expression A - B
  # i.e. our DE analyses examine A wrt B - CHANGED

  ##################################
  # TODO: testing
  v2 <- as.character(comparison.pairs$V2)
  for (i in 1:length(v2)) {
      if (!is.na(as.numeric(substr(v2[i],1,1)))) {
          v2[i] <- paste0('X', v2[i])
      }
  }
  comparison.pairs$V2 <- v2
  ##################################
  comparison.exprs = paste(comparison.pairs$V1, comparison.pairs$V2, sep='-')
  #  print(comparison.exprs)
  contra = makeContrasts(contrasts=comparison.exprs, levels=design)
  colnames(contra) = comparison.names
  return (contra)
}

run_glm_batch_contrasts <- function(counts, design, contrasts, out.dir) {
  dge = DGEList(counts, remove.zeros=T)
  dge = calcNormFactors(dge, method='TMM')

  # check to see if and groups have only one replicate -- if so, perform modified modeling
  if (any(colSums(design)==1)) {
      dge = estimateGLMCommonDisp(dge, method = 'deviance', robust = TRUE, subset = NULL)
      dge = estimateGLMTrendedDisp(dge, method = 'auto')
      dge = estimateGLMTagwiseDisp(dge)
  } else {
      dge = estimateGLMCommonDisp(dge, design)
      dge = estimateGLMTrendedDisp(dge, design)
      dge = estimateGLMTagwiseDisp(dge, design)
  }
  fit = glmFit(dge, design)
  for (ii in 1:ncol(contrasts)) {
    comparison.name = colnames(contrasts)[ii]
    comparison.coeffs = as.numeric(contrasts[,ii])
    print(comparison.coeffs)
    lrt = glmLRT(fit, contrast=comparison.coeffs)
    ranked = topTags(lrt, n=NULL, sort.by = 'PValue')

    # TODO: add Lingfei suggestions here?
    # for each gene in ranked, add logCPM for each

    out.file = file.path(out.dir, paste(comparison.name, 'edgeR', 'tsv', sep='.'))
    write.table(ranked, file=out.file, sep='\t', quote=F, row.names=T)
  }
}

get_samples <- function(samples.described.file) {
  samples = read.delim(samples.described.file, row.names=2, header=FALSE)
  return(samples)
}

get_comparisons <- function(samples.compared.file) {
  comparison.names = make.names(readLines(samples.compared.file))
  comparison.pairs = gsub('_vs_', ' ', comparison.names)
  comparison.pairs = read.table(textConnection(comparison.pairs))
  return (comparison.pairs)
}

run_classic_batch_contrasts <- function(counts, samples, comparisons, out.dir) {
  for (ii in 1:nrow(comparisons)) {
    a.class = as.character(comparisons[ii,1])
    b.class = as.character(comparisons[ii,2])
    a.cols = rownames(samples)[which(samples$V1 == a.class)]
    b.cols = rownames(samples)[which(samples$V1 == b.class)]
    comparison.name = paste(a.class, b.class, sep='_vs_')
    out.file = file.path(out.dir, paste(comparison.name, 'edgeR', 'tsv', sep='.'))
    run_edger_classic(counts, a.cols, b.cols, out.file=out.file)
  }
}

combine_edgeR_results <- function(out_dir) {
    # get list of files
    individual_files <- file.path(out_dir, list.files(out_dir))
    individual_files <- individual_files[str_detect(individual_files, '\\.edgeR\\.tsv')]
    individual_files <- individual_files[basename(individual_files) != 'combined_edgeR.tsv']

    individual_file_list <- vector('list', length = length(individual_files))
    for (i in 1:length(individual_files)) {
        individual_file_list[[i]] <- fread(individual_files[i], sep = '\t')
        comparison_name <- str_replace(basename(individual_files[i]), '.edgeR.tsv', '')
        individual_file_list[[i]][, comparison := comparison_name]
    }
    combined_edgeR <- rbindlist(individual_file_list)
    setnames(combined_edgeR, 'V1', 'gene')
    combined_edgeR[, FDR := p.adjust(PValue, method = 'BH')]
    combined_edgeR <- setkey(combined_edgeR, by = 'FDR')
    fwrite(combined_edgeR, file.path(out_dir, 'combined_edgeR.tsv'), sep = '\t')
}

# main
main <- function() {
    args = commandArgs(trailingOnly=TRUE)
    counts_file <- args[1]
    samples_file <- args[2]
    comparisons_file <- args[3]
    out_dir <- args[4]
    dir.create(out_dir, showWarnings = FALSE)

    counts = read.csv(counts_file, row.names=1)
    # remove ensembl_id
    counts <- subset(counts, select = -c(ensembl_id))

    design = get_simple_design(samples_file)
    rownames(design) <- make.names(str_trim(rownames(design)))

    if (!setequal(names(counts), rownames(design))) {
        names(counts) <- str_replace_all(names(counts), '_S[0-9]+', '')
    }

    counts.new = counts[, make.names(rownames(design))]

    contrasts = get_simple_contrasts(comparisons_file, design)

    # filter counts to only keep genes with some real activity across samples
    # (in this case, keep the top 8000 genes in terms of cross-sample variance of CPM)
    #counts = filter_cpm_rank(counts, keep.genes=4000, method='var')

    # remove genes with zero counts across all samples
    counts.new = filter_zeros(counts.new)
    counts.new <- na.omit(counts.new)

    # run DE analysis
    run_glm_batch_contrasts(counts.new, design, contrasts, out_dir)

    # combine all edgeR results
    combine_edgeR_results(out_dir)

}

main()
