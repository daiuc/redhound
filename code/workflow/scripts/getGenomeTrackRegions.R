

#---------------------------argparse---------------------------------

# argparse requires python in the same environment as R
suppressPackageStartupMessages(library(argparse, quietly = T))

parser <- ArgumentParser()

parser$add_argument("-I", "--colocInput", dest = "colocInput", type = "character",
                    required = T, help = "Annotated coloc summary file.")

parser$add_argument("-R", "--regions", dest = "regions", type = "character",
                    required = T, help = "locusCompare ploted phenotype ids")

parser$add_argument("-O", "--fout", dest = "fout", type = "character",
                    required = T, help = "out file")


args <- parser$parse_args()

coloc.f <- args$colocInput
regions.f <- args$regions
fout <- args$fout

#----------------------library-----------------------------

suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(purrr))
suppressMessages(library(GenomicRanges))


if (interactive()) {
    coloc.f <- "results/chRNA/YRI/coloc/sQTL/coloc_summary_annotated.tsv"
    regions.f <- "resources/230220-edQTL-sQTL-LocusComp-Handpicked-ids.txt"
    fout <- "results/chRNA/YRI/coloc/sQTL/coloc_for_pygenometrack.tsv"
}


#----------------------functions-----------------------------

concatLabels <- function(data, byCol, oldLabel, newLabel) {
    #' @param data    : data.table
    #' @param byCol   : group by column, typically key
    #' @param oldLabel: label to concat
    #' @param newLabel: label after concat

    env = list(byCol = byCol, oldLabel = oldLabel, newLabel = newLabel)
    env = lapply(env, as.name)
    expr1 = substitute(paste0(unique(oldLabel), collapse = ","), env)
    expr2 = substitute(byCol, env)

    data = eval(substitute(data[, .(new = expr1), by = expr2]))
    setnames(data, "new", newLabel)
    unique(data)
}


#----------------------main-----------------------------

# annotated coloc results
coloc <- fread(coloc.f, colClasses = c(topsnpsin = "character",
    edsitein = "character", edsiteinNC = "character"))
coloc[, id := paste(pheno1, pheno2, sep="|")]

# hand picked regions
regions <- readLines(regions.f)
regions <- str_replace_all(regions, "[\\-\\−]", "-") %>%
    str_replace_all("\\ ", "")

if (!all(regions %in% coloc$id)) {
  stop(glue:glue("Error! Not all keys in {regions.f} are in {coloc.f}!"))
}

cat(glue::glue("\n\nAre all regions in coloc summary? {all(regions %in% coloc$id)}\n\n"))

# filter coloc to select regions only
coloc <- coloc[id %in% regions]

# fix multiple annotation per id
dt1 <- concatLabels(coloc, "id", "genename", "genename")
dt2 <- concatLabels(coloc, "id", "genefeature", "genefeature")
dt3 <- concatLabels(coloc, "id", "genetype", "genetype")

coloc <-  coloc[, -c("genename", "genefeature", "genetype")] %>% unique
coloc <- coloc[dt1, on="id", nomatch=NULL] %>%
    .[dt2, on="id", nomatch=NULL] %>%
    .[dt3, on="id", nomatch=NULL]

# plotted regions
gr <- coloc[, .(p1window, p2window)] %>% map(~as(.x, "GRanges"))

gr <- map(1:length(gr[[1]]),
    ~GenomicRanges::reduce(c(gr$p1window[.x], gr$p2window[.x]))) %>%
    purrr::reduce(c)
mcols(gr) <- coloc


cat("\n\nWrite coloc regions for pygenometracks.\n\n")
fwrite(as.data.table(gr), fout, sep = "\t")
