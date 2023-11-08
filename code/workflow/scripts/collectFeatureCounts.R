#!/usr/bin/env Rscript

#'---------------------------------------------------------------
#' @version     : v1
#' @Date        : 2023-02-22
#' @author      : Chao Dai
#' @Description : collect featureCounts output into a single count table
#'---------------------------------------------------------------


#---------------------------argparse---------------------------------

# argparse requires python in the same environment as R
suppressPackageStartupMessages(library(argparse, quietly = T))

parser <- ArgumentParser()

parser$add_argument("-C", "--countFiles", dest = "countFiles", type = "character",
    required = T, help = "count files from featureCounts. comma separated list.")

parser$add_argument("-N", "--nskip", dest = "nskip", type = "integer",
    required = F, default = 2, help = "nrows to skip. default 2.")

parser$add_argument("-O", "--fout", dest = "fout", type = "character",
    required = T, help = "out file, tsv file with no header")

args       <- parser$parse_args()
countFiles <- args$countFiles
nSkip      <- args$nskip
fout       <- args$fout

if (interactive()) {
    countFiles <- c("results/mRNA/YRI/featureCounts/NA18486.counts",
                    'results/mRNA/YRI/featureCounts/NA18487.counts',
                    'results/mRNA/YRI/featureCounts/NA19119.counts')
    nSkip <- 2
}


suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(stringr))


countFiles <- str_split(countFiles, "[,\\ ]") %>% unlist
samples <- map_chr(countFiles, ~str_extract(.x, 'NA\\d+'))
names(countFiles) <- samples

cat(glue::glue("\n\nCollecting counts from {length(samples)} sample:\n\n"))
print(samples)

# read all count tables
counts <- map(countFiles, ~fread(.x, header = F, skip = nSkip))
geneCols <- map(counts, ~.x$V1) # gene_id column for each table
valCols <- map(counts, ~.x$V7) # count value columns

geneLen <- counts[[1]][, .(gene_id = V1, gene_length = V6)]

if (all(purrr::reduce(geneCols, intersect) == geneCols[[1]])) {
    cat("\n\nCombining gene expression counts into one count table.\n\n")
    dt <- as.data.table(valCols)
    dt <- dt[, c(list(gene_id = geneCols[[1]]), .SD)]

    # for expressed gene list, require on average at least half samples have
    # at least 1 read
    row_sums <- dt[, -c("gene_id")] %>% rowSums
    keep_rows <- row_sums > (ncol(dt) - 1) * .5
    keep_genes <- dt[keep_rows, gene_id]

    gout <- file.path(dirname(fout), "expressed_genes.tsv")
    cat(glue::glue("\nWrite list of expressed genes to {gout}"))
    writeLines(keep_genes, gout)

    cat(glue::glue("\n\nWrite to file: {fout}..."))
    fwrite(dt, fout, col.names = T, sep = "\t")

} else {
    stop("Error! featureCounts output files must have consistent gene_id across files!")
}
