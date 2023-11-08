#!/usr/bin/env Rscript

#---------------------------------------------------------------
#
#' @author      : Chao Dai
#' @Description : Use this script to run coloc on molecular traits
#' 
#
#---------------------------------------------------------------



#---------------------------argparse---------------------------------

# argparse requires python in the same environment as R
suppressPackageStartupMessages(library(argparse, quietly = T))

parser = ArgumentParser()

parser$add_argument("--perm1", dest = "perm1", type = "character",
    required = T, help = "perm pass result for mol trait 1")

parser$add_argument("--perm1cols", dest = "perm1cols", type = "character",
    required = T, help = "colnames in the trait 1 permuation result file, \
    separted by coma. e.g. pid,pchr,pstart")

parser$add_argument("--nom1", dest = "nom1", type = "character",
    required = T, help = "nominal pass result for mol trait 1. May be a file \
    or prefix")

parser$add_argument("--nom1_suf", dest = "nom1_suf", type = "character",
    required = F, help = "If --nom1 is prefix, supply suffix here. Structure \
    is `{prefix}/{chrom}/{suffix}`")

parser$add_argument("--nom1cols", dest = "nom1cols", type = "character",
    required = T, help = "colnames in the trait 1 nominal pass result file, \
    separted by coma. e.g. pid,pchr,pstart")

parser$add_argument("--perm2", dest = "perm2", type = "character",
    required = T, help = "perm pass result for mol trait 2")

parser$add_argument("--perm2cols", dest = "perm2cols", type = "character",
    required = T, help = "colnames in the trait 2 permuation result file, \
    separted by coma. e.g. pid,pchr,pstart. Make sure")

parser$add_argument("--nom2", dest = "nom2", type = "character",
    required = T, help = "nominal pass result for mol trait 2. May be a file \
    or prefix")

parser$add_argument("--nom2_suf", dest = "nom2_suf", type = "character",
    required = F, help = "If --nom2 is prefix, supply suffix here. Structure is\
     `{prefix}/{chrom}/{suffix}`")

parser$add_argument("--nom2cols", dest = "nom2cols", type = "character",
    required = T, help = "colnames in the trait 2 nominal pass result file, \
    separted by coma. e.g. pid,pchr,pstart")

parser$add_argument("-m", "-minOverlap", dest = "minOverlap", default = 50000,
    type = "integer", help = "minimum bases of overlap between trait1 & trait2, default 50kb.")

parser$add_argument("--extend1", dest = "extend1", default = 100000,
    type = "integer", help = "extend trait 1 by extend1 on both end. Default 100kb.")

parser$add_argument("--extend2", dest = "extend2", default = 0,
    type = "integer", help = "extend trait 2 by extend2 on both end. Default 0.")

parser$add_argument("-O", "--fout", dest = "fout", type = "character", 
    required = T, help = "out file")

parser$add_argument("-c", "--core", dest = "core", type = "integer",
    default = 1, help = "Number of core to use, default 1.")

args = parser$parse_args()



#---------------------------library------------------------------

suppressMessages(library(stringr))
suppressMessages(library(purrr))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(glue))
suppressMessages(library(coloc))
suppressMessages(library(furrr))
suppressPackageStartupMessages((library(tictoc, quietly = T)))




#---------------------------functions------------------------------

readNominal = function(nominal.f, nominal.suf = "Nominal.txt.gz", gr, nom_cols) {
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

  c = seqnames(gr) %>% as.vector
  s = start(gr)
  e = end(gr)

  if (file.info(nominal.f)$isdir) {
    # fill in full file path using chrom from gr, and prefix, suffix
    f = glue(nominal.f, "/", {c}, "/", {nominal.suf})
  } else {
    f = nominal.f
  }
  if (!file.exists(f)) stop(glue("File {f} does not exist!"))

  CMD = glue("tabix {f} {c}:{s}-{e}")
  dt = fread(cmd = CMD, col.names = nom_cols) # read nominal snps in  region
  if (nrow(dt) > 0) {
    dt = dt[pid %in% mcols(gr)[['pid']]] # only want snps for given pheno
  }
  return(dt)
}


getDataset = function(gr, nom.file, nom_cols) {
  #' @description     : prepare qtl nominal pass dataset for coloc test
  #'
  #' @param gr        : GRange object, of the coloc window for the phenotype to
  #'                    be tested
  #' @param nom.file  : character, path (or prefix) to the nominal pass result
  #' @param nom_cols  : character, column names for nominal pass result. Make
  #'                    sure this
  #'                    is consistent across dataset
  #'
  #' @return          : data.table, dataset with necessary columns for coloc
  #'                    dataset. Note since this dataset does not ensure all
  #'                    snps are shared and sorted identically with the other
  #'                    coloc dataset

  pid = mcols(gr)[['pid']] # one phenotype id

  # read nominal pass output for the phenotype, use tabix
  ds = readNominal(nominal.f = nom.file, gr = gr, nom_cols = nom_cols)

  if (nrow(ds) < 10) {
    cat(glue("\n\nNot enough SNPs for {pid}. Return NULL."))
    return(NULL)
  } else {
    # compute z and varbeta
    ds[, z := qnorm(pval.nom / 2, lower.tail = F)]
    ds[, varbeta := (beta / z)^2]
    ds = ds[!is.na(varbeta)] # remove NA

    # some indels have multiple rows with same snp_id, keep 1
    ds = ds[order(pval.nom), .SD, by = c("pid", "snp_id")] %>%
        unique(by = c("pid", "snp_id"))

    return(ds)
  }

}


hush = function(code) {
    # shut up coloc. credit: stackoverflow
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
}






#---------------------------setup---------------------------------


perm1.f = args$perm1
nom1.f = args$nom1
if (!is.null(args$nom1_suf)) nom1_suf = args$nom1_suf
perm1.cols = str_split(args$perm1cols, ",")  %>% unlist
nom1.cols = str_split(args$nom1cols, ",")  %>% unlist

perm2.f = args$perm2
nom2.f = args$nom2
if (!is.null(args$nom2_suf)) nom2_suf = args$nom2_suf
perm2.cols = str_split(args$perm2cols, ",")  %>% unlist
nom2.cols = str_split(args$nom2cols, ",")  %>% unlist

minOverlap = args$minOverlap
extend1 = args$extend1
extend2 = args$extend2

# parallel setup future_map
ncore = args$core
ncore = min(ncore, availableCores())
cat(glue("\n\nRunning with {ncore} cores.\n\n"))
plan(multisession, workers = ncore)

# output
fout = args$fout

cat(glue("\n\nPermuation pass for trait 1:\n    {perm1.f}\n"))
cat(glue("\n\nNominal pass for trait 1:\n    {nom1.f}\n"))
cat(glue("\n\nPermuation pass for trait 2:\n    {perm2.f}\n"))
cat(glue("\n\nNominal pass for trait 2:\n    {nom2.f}\n\n"))

cat(glue("\n\nColoc tests require trait 1 and trait 2 has min {minOverlap} bp \
          overlap.\n\n"))


#---------------------------Load---------------------------------


if (!("pval_beta_adj" %in% perm1.cols) | !("pval_beta_adj" %in% perm2.cols)) {
    stop("Stop! Make sure'pval_beta_adj' is the column name for beta adjusted pvalue")
}

cat("\n\nLoading permutation pass results...")
perm1 = fread(perm1.f, col.names = perm1.cols)
perm2 = fread(perm2.f, col.names = perm2.cols)

cat("\n\nRemoving qtls with beta adjusted pval >= 0.01")
perm1 = perm1[pval_beta_adj < 0.01]
perm2 = perm2[pval_beta_adj < 0.01]

cat(glue("\n\nqtl1 has {nrow(perm1)} phenotypes."))
cat(glue("\nqtl2 has {nrow(perm2)} phenotypes."))


#---------------------------Find Overlaps--------------------------

perm1.gr = makeGRangesFromDataFrame(
  perm1, keep.extra.columns = T, ignore.strand = F,
  seqnames.field = "pchr",
  start.field = "pstart",
  end.field = "pend",
  strand.field = "pstrand"
)
cat(glue("\n\nExtend {extend1} both end for trait1's phenotypes."))
perm1.gr = perm1.gr + extend1

perm2.gr = makeGRangesFromDataFrame(
  perm2, keep.extra.columns = T, ignore.strand = F,
  seqnames.field = "pchr",
  start.field = "pstart",
  end.field = "pend",
  strand.field = "pstrand"
)

cat(glue("\n\nExtend {extend2} both end for trait2's phenotypes."))
perm2.gr = perm2.gr + extend2

# overlaps
olaps = findOverlaps(perm1.gr, perm2.gr, minoverlap = minOverlap)
glue("\n\nFound {length(olaps)} pairs of overlapping regions between trait 1 and trait 2.")


#---------------------------coloc-------------------------------

n_max = length(olaps)
# n_max = 30 # test

if (n_max < length(olaps)) cat(glue("\n\nTesting! Running {n_max} tests.\n\n"))

tic()

res.l = future_map(
  .x = 1:n_max,
  .f = function(x) {

    if (x %% 500 == 0) {
        t = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        cat(glue("\n\n{t}...Number of tests run: ~{x}"))
    }

    o = olaps[x] # a single findOverlaps object
    p1.gr = perm1.gr[queryHits(o)]
    p2.gr = perm2.gr[subjectHits(o)]
    snps1 = getDataset(p1.gr, nom1.f, nom1.cols)
    snps2 = getDataset(p2.gr, nom2.f, nom2.cols)

    if (!is.null(snps1) & !is.null(snps2)) {

        # subset only shared snps and ensure same order
        ss = intersect(snps1$snp_id, snps2$snp_id)

        # at least 20+ shared snps
        if (length(ss) > 20) {

            snps1 = snps1[snp_id %in% ss][order(gstart)]
            snps2 = snps2[snp_id %in% ss][order(gstart)]

            ds1 = list(
                beta = snps1$beta,
                varbeta = snps1$varbeta,
                type = "quant",
                sdY = 1,
                snp = snps1$snp_id,
                position = snps1$gstart
            )

            ds2 = list(
                beta = snps2$beta,
                varbeta = snps2$varbeta,
                type = "quant",
                sdY = 1,
                snp = snps1$snp_id,
                position = snps1$gstart
            )

            # coloc test
            res = suppressWarnings(hush(coloc.abf(ds1, ds2)))
            res = res$summary %>% t %>% as.data.table
            res[, `:=`(
                pheno1 = mcols(p1.gr)[["pid"]],
                pheno2 = mcols(p2.gr)[["pid"]],
                p1window = as.character(p1.gr),
                p2window = as.character(p2.gr),
                olapwindow = as.character(intersect(p1.gr, p2.gr)),
                p1topsnp = mcols(p1.gr)[["best_gid"]],
                p2topsnp = mcols(p2.gr)[["best_gid"]]
                )]
            setorder(res, -`PP.H4.abf`)
            return(res)
        }
    }
  }
)

res.l = rbindlist(res.l)
setorder(res.l, -`PP.H4.abf`)

cat("\n\n")
toc()
cat(glue("\n\nTotal tests run: ~{n_max}"))

cat(glue("\n\nWriting result to {fout}..."))
fwrite(res.l, fout, sep = "\t")
cat("Done.")


