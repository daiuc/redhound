#!/usr/bin/env Rscript

#'---------------------------------------------------------------
#' @version     : v2
#' @Date        : 2023-02-17
#' @author      : Chao Dai
#' @Description : Plot locusCompare plots for coloc results
#' @note        : locuscomparer requires internet connection!
#'---------------------------------------------------------------


#---------------------------argparse---------------------------------

# argparse requires python in the same environment as R
suppressPackageStartupMessages(library(argparse, quietly = T))

parser <- ArgumentParser()

parser$add_argument("--colocFile", dest = "colocFile", type = "character",
                    required = T, help = "Annotated coloc summary file.")

parser$add_argument("--nom1", dest = "nom1", type = "character",
                    required = T, help = "Nominal file or prefix")

parser$add_argument("--nom1cols", dest = "nom1cols", type = "character",
                    required = T, help = "Nominal pass column names, comma delimited.")

parser$add_argument("--nom2", dest = "nom2", type = "character",
                    required = T, help = "Nominal file or prefix")

parser$add_argument("--nom2cols", dest = "nom2cols", type = "character",
                    required = T, help = "Nominal pass column names. comma delimited.")

parser$add_argument("--minH4PP", dest = "minH4PP", type = "double",
                    default = 0.7, help = "Minimum H4 PP, default 0.7.")

parser$add_argument("-O", "--fout", dest = "fout", type = "character",
                    required = T, help = "out file")

parser$add_argument("-c", "--core", dest = "core", type = "integer",
                    default = 1, help = "Number of core to use, default 1.")

args <- parser$parse_args()

coloc.f <- args$colocFile
nom1.f <- args$nom1
nom2.f <- args$nom2
nom1cols <- args$nom1cols
nom2cols <- args$nom2cols
pp4thd <- args$minH4PP
fout <- args$fout







#-------------------------functions---------------------------------

readNominal <- function(nominal.f, nominal.suf = "Nominal.txt.gz", gr, nom_cols) {
  #' @description  : read nominal pass results into a data.table
  #'
  #' @nominal.f    : character, tabix indexed nominal pass result file prefix
  #' @nominal.suf  : character, suffice for the nominal file, if nominal.f is
  #'                 given as prefix, default "Nominal.txt.gz"
  #' @gr           : GRange object, of the queried region
  #' @pid_col      : character, column name for phenotype id
  #' @nom_cols     : vector, column nmaes for nominal pass file
  #'
  #' @return       : data.table, nominal pass results for the queried region
  #'                 for the phenotype

  c <- seqnames(gr) %>% as.vector
  s <- start(gr)
  e <- end(gr)

  if (file.info(nominal.f)$isdir) {
    # fill in full file path using chrom from gr, and prefix, suffix
    f <- glue(nominal.f, "/", {c}, "/", {nominal.suf})
  } else {
    f <- nominal.f
  }
  if (!file.exists(f)) stop(glue("File {f} does not exist!"))

  CMD <- glue("tabix {f} {c}:{s}-{e}")
  dt <- fread(cmd = CMD, col.names = nom_cols) # read nominal snps in  region
  if (nrow(dt) > 0) {
    dt <- dt[pid %in% mcols(gr)[['pid']]] # only want snps for given pheno
  }

  return(dt)
}


getLocusCompareInputs <- function(colocid, colocdt, snpdb,
                                  nom1file, nom1cols,
                                  nom2file, nom2cols) {
  #' @param colocid : colocalization id (from dt)
  #' @param colocdt : colocalization summary dt
  #' @param snpdb   : dbSNP object from import
  #' @param nom1file: nominal path output file for trait 1
  #' @param nom1cols: column names for nom1file
  #' @param nom2file: nominal path output file for trait 2
  #' @param nom2cols: column names for nom2file
  #'
  #' @return : write 2 temp file, one each for trait 1 and trait 2.
  #'           Out file has two columns: rsid, pval. Ready for plottting
  #'           locusCompare plots.
  #'
  #'

  # select coloc test based on id
  c1 = colocdt[id == colocid][1] # in case of multiple records

  # merge pheno1 and pheno2 windows
  gr1 <- reduce(c(as(c1$p1window, "GRanges"), as(c1$p2window, "GRanges")))
  gr2 <- gr1 # same GRanges, but differenet pheno id
  mcols(gr1)$pid <- c1$pheno1
  mcols(gr2)$pid <- c1$pheno2

  # read in all nominal pass SNPS in the region for each trait
  snps1 <- readNominal(nominal.f = nom1file, gr = gr1, nom_cols = nom1cols)
  snps2 <- readNominal(nominal.f = nom2file, gr = gr2, nom_cols = nom2cols)

  # grap all snps in dbsnp for the merged region
  gr <- gr1
  seqlevelsStyle(gr) <- "NCBI"
  rsids <- snpsByOverlaps(x = snpdb, ranges = gr)

  # snps from trait 1 and trait 2 as GRanges
  snps1 <- makeGRangesFromDataFrame(
    snps1, keep.extra.columns = T, ignore.strand = T,
    seqnames.field = "gchr", start.field = "gstart",
    end.field = "gend"
  )
  seqlevelsStyle(snps1) <- "NCBI"

  snps2 <- makeGRangesFromDataFrame(
    snps2, keep.extra.columns = T, ignore.strand = T,
    seqnames.field = "gchr", start.field = "gstart",
    end.field = "gend"
  )
  seqlevelsStyle(snps2) <- "NCBI"

  # annotate trait 1 snps with dbsnp rsid
  o1 = findOverlaps(snps1, rsids)
  snps1= snps1[queryHits(o1)]
  mcols(snps1) = cbind(mcols(snps1), mcols(rsids[subjectHits(o1)]))
  seqlevelsStyle(snps1) <- "UCSC"

  # annotate trait 2 snps with dbsnp rsid
  o2 = findOverlaps(snps2, rsids)
  snps2= snps2[queryHits(o2)]
  mcols(snps2) = cbind(mcols(snps2), mcols(rsids[subjectHits(o2)]))
  seqlevelsStyle(snps2) <- "UCSC"

  # to data.table
  snps1 <- as.data.table(snps1)
  setnames(snps1, c("RefSNP_id", "pval.nom"), c("rsid", "pval"))
  snps2 = as.data.table(snps2)
  setnames(snps2, c("RefSNP_id", "pval.nom"), c("rsid", "pval"))

  # remove duplicates due to one snp having multiple rsid
  snps1 <- unique(snps1, by = "rsid")
  snps2 <- unique(snps2, by = "rsid")

  return(list(trait1 = snps1, trait2 = snps2))
}


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
  return(unique(data))
}



#---------------------------library---------------------------------



suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(glue))
suppressMessages(library(GenomicRanges))
suppressMessages(library(SNPlocs.Hsapiens.dbSNP151.GRCh38))
library(locuscomparer)


#---------------------------debug---------------------------------

if (interactive()) {
  print("Running in interactive mode!!\n")
  coloc.f <- "results/chRNA/YRI/coloc/eQTL/coloc_summary_annotated.tsv"

  nom1.f <- "results/chRNA/YRI/edQTL_FixedPCs/ranknorm/cis_1000000/"
  nom2.f <- "resources/bjf79-qtltools/chRNA.Splicing/NominalPassForColoc.txt.tabix.gz"

  nom1cols <- str_split('id,pchr,pstart,pend,pstrand,n_snps,dist,snp_id,gchr,gstart,gend,pval.nom,pval.r2,beta,topSNP', ",") %>% unlist
  nom2cols <- str_split('pid,pchr,pstart,pend,pstrand,n_snps,dist,snp_id,gchr,gstart,gend,pval.nom,pval.r2,beta,beta_se,topSNP', ",") %>% unlist

  pp4thd <- .7
  fout  <- "test/test.eQTL.plotLocus.json"
}



#---------------------------load---------------------------------

cat("\n\nLoad annotated coloc summary...\n\n")
coloc <- fread(coloc.f)  %>% .[PP.H4.abf > pp4thd]
setorder(coloc, -PP.H4.abf)
coloc[, id := paste(pheno1, pheno2, sep = "|")]

nom1cols <- str_split(nom1cols, ",") %>% unlist
nom2cols <- str_split(nom2cols, ",") %>% unlist

# dbSNP151 db
cat("\n\nGet dbSNP151 database.\n\n")
dbSNP <- SNPlocs.Hsapiens.dbSNP151.GRCh38


#---------------------------main---------------------------------

# ----------- Ensure 1 id has 1 record ---------------
coloc2 <- map(c("genename", "genefeature", "genetype", "repeat", "db"),
  ~concatLabels(coloc, "id", .x, .x)) %>%
  purrr::reduce(., dplyr::inner_join, by = c("id"))

cols1 <- c("nsnps", "PP.H4.abf", "pheno1", "pheno2", "p1window", "p2window",
          "olapwindow", "p1topsnp", "p2topsnp", "topsnpsin", "edsitein",
          "edsiteinNC", "id")
coloc <- coloc[, ..cols1] %>% unique
coloc = coloc[coloc2, on = "id", nomatch = NULL]

# ----------- Get LocusCompare input data ---------------

# test >>>
# ids <- coloc$id[1:5]
# names(ids) <- ids
# <<<

ids <- unique(coloc$id)
names(ids) <- ids

cat(glue("\n\nPrepare rsid and pval for {length(ids)} coloc tests...\n\n"))
tictoc::tic()
dts <- base::lapply(ids, getLocusCompareInputs,
              colocdt = coloc, snpdb = dbSNP,nom1file = nom1.f,
              nom1cols = nom1cols, nom2file = nom2.f, nom2cols = nom2cols)

tictoc::toc()

# convert data.table to json and save to file
print("Convert and write locusCompare input data into json file.")
dts <- lapply(dts, function(ll) {
  l <- map(ll, ~ jsonlite::toJSON(.x, "column"))
  return(l)
})

jsonlite::write_json(dts, fout)




