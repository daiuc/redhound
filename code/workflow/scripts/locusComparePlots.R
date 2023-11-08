#!/usr/bin/env Rscript

#'---------------------------------------------------------------
#' @version     : v1
#' @Date        : 2023-3-3
#' @author      : Chao Dai
#' @Description : Plot locusCompare plots for coloc results
#' @note        : locuscomparer requires internet connection!
#'---------------------------------------------------------------




#---------------------------argparse---------------------------------

# argparse requires python in the same environment as R
suppressPackageStartupMessages(library(argparse, quietly = T))

parser <- ArgumentParser(description = "plot locuscompare plots.")

parser$add_argument("--jsonFile", dest = "jsonFile", type = "character",
                    required = T, help = "json file for locusCompare data")

parser$add_argument("--colocFile", dest = "colocFile", type = "character",
                    required = T, help = "Annotated coloc summary file.")

parser$add_argument("--perm1", dest = "perm1", type = "character",
                    required = T, help = "permutation pass for trait 1 (edQTL)")

parser$add_argument("--perm1cols", dest = "perm1cols", type = "character",
                    required = T, help = "permutation pass columns for trait 1. comma delimited.")

parser$add_argument("--perm2", dest = "perm2", type = "character",
                    required = T, help = "permutation pass for trait 2 (eQTL, sQTL, ncRNA QTL)")

parser$add_argument("--perm2cols", dest = "perm2cols", type = "character",
                    required = T, help = "permutation pass columns for trait 2. comma delimited.")

parser$add_argument("-O", "--outFile", dest = "outFile", type = "character",
                    required = T, help = "output.")


args <- parser$parse_args()

json.file = args$jsonFile
colocFile = args$colocFile
perm1File = args$perm1
perm2File = args$perm2
perm1cols = args$perm1cols
perm2cols = args$perm2cols
fout = args$outFile


if (interactive()) {
    print("Interactive mode.")
    json.file = 'results/chRNA/YRI/coloc/eQTL/locusCompareData.json'
    colocFile = 'results/chRNA/YRI/coloc/eQTL/coloc_summary_annotated.tsv'
    perm1File = 'results/chRNA/YRI/edQTL_FixedPCs/ranknorm/cis_10000/permutations_pass_all_chr.AddedQvalue.txt.gz'
    perm2File = 'resources/bjf79-qtltools/chRNA.Expression.Splicing/PermutationPassForColoc.txt.gz'
    perm1cols = 'pid,pchr,pstart,pend,pstrand,n_snps,best_dist,best_gid,best_gchr,best_gstart,best_gend,dof_true,dof_est,beta_ml1,beta_ml2,pval_nom,pval_r2,slope,pval_emp,pval_beta_adj,q'
    perm2cols = 'pid,pchr,pstart,pend,pstrand,n_snps,best_dist,best_gid,best_gchr,best_gstart,best_gend,dof_true,dof_est,beta_ml1,beta_ml2,pval_nom,pval_r2,slope,slope_se,pval_emp,pval_beta_adj'
    fout = 'results/chRNA/YRI/coloc/eQTL/plot-test.pdf'
}



# ---------------- Setup ---------------------

suppressMessages(library(jsonlite))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(glue))
suppressMessages(library(locuscomparer))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))



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




# ---------------- coloc summary ---------------------

coloc = fread(colocFile)
coloc[, id := paste(pheno1, pheno2, sep = "|")]

# ----------- Ensure 1 id has 1 record ---------------
coloc2 <- map(c("genename", "genefeature", "genetype", "repeat", "db"),
              ~concatLabels(coloc, "id", .x, .x)
             ) %>%
            reduce(., dplyr::inner_join, by = c("id"))

cols1 <- c("nsnps", "PP.H4.abf", "pheno1", "pheno2", "p1window", "p2window",
          "olapwindow", "p1topsnp", "p2topsnp", "topsnpsin", "edsitein",
          "edsiteinNC", "id")
coloc <- coloc[, ..cols1] %>% unique
coloc = coloc[coloc2, on = "id", nomatch = NULL]


# ---------------- read in json for plotting ---------------------
print("Read in json data file for locusCompare SNP data")
plotdata = read_json(json.file, simplifyVector = T)
plotdata = lapply(
    plotdata,
    function(ll) {
        ll = map(ll, ~ fromJSON(.x)  %>% as.data.table)
    }
)

ids <- names(plotdata) # coloc test ids, trait1_pid|trait2_pid
names(ids) <- ids

# trait 1 pid, trait2 pid for each test id. used to get betas next
pids <- lapply(
  ids,
  function(x) {
    s <- str_split(x, "\\|") %>% unlist
    return(list(pheno1 = s[1], pheno2 = s[2]))
  }
)


# ---------------- get betas from perm pass ---------------------

perm1cols <- str_split(perm1cols, ",")  %>% unlist
perm2cols <- str_split(perm2cols, ",")  %>% unlist

perm1 <- fread(perm1File, col.names = perm1cols)
perm2 <- fread(perm2File, col.names = perm2cols)

betas <- lapply(
  pids,
  function(ll) {
    bt <- map2(.x = ll, .y = list(perm1, perm2),
               ~.y[pid == .x, slope])
    bt <- glue("trait1-beta: {bt$pheno1}, trait2-beta: {bt$pheno2}")
    return(bt)
  }
)

rm(perm1, perm2)
gc()


# ---------------- make plots ---------------------

# maxN  <- 10 # <<< for testing
maxN <- length(plotdata)

print(glue("Generating {length(plotdata[1:maxN])} LocusCompare plots."))

# population must be one of the 5 popuations from 1000 Genomes: 
# 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'

ggs <- map(plotdata[1:maxN],
  ~locuscompare(.x$trait1[, .(rsid, pval)],
                .x$trait2[, .(rsid, pval)],
                title1 = "edQTL",
                title2 = "eQTL",
                population = "AFR",
                genome = "hg38")
  )

# add titles
print("Adding details to plot titles.")
plot_titles = map_chr(ids[1:maxN],
  ~coloc[id == .x] %>%
    .[1, .(t = paste("trait1:", pheno1, "trait2:", pheno2,
                     "\nBothTopSNPsOlap?(bits):", topsnpsin,
                     "A2IInOlap?(bits):", edsitein, "A2IInNCRNA:", edsiteinNC,
                     "\nA2IF:", genefeature, "A2IInGene:", genename,
                     "A2IIsRepeat:", `repeat`, "A2IInDB:", db))
      ] %>% .$t
)

# add in betas to titles
plot_titles <- paste(plot_titles, betas[1:maxN], sep = '\n')


# ---------------- make plots ---------------------

cat("\n\nCombining plots to PDF, 1 plot per page.\n\n")
cat(glue("\n\nSaving plots to {fout}\n"))

# pdf(fout, width=8, height=5)
ggs.togo = marrangeGrob(
  grobs = ggs,
  nrow=1,
  ncol=1,
  top = quote(plot_titles[g])
)

ggsave(fout, ggs.togo, device="pdf", width=8, height=7)

cat(glue("\n\n{Sys.time()} Done.\n"))


