#!/usr/bin/env Rscript

#---------------------------------------------------------------
#
#' @author      : Chao Dai
#' @Description : Annotate coloc results
#' 
#
#---------------------------------------------------------------


#---------------------------argparse---------------------------------



# argparse requires python in the same environment as R
suppressPackageStartupMessages(library(argparse, quietly = T))

parser = ArgumentParser()

parser$add_argument("-I", "--coloc", dest = "coloc", type = "character",
    required = T, help = "colocalization result")

parser$add_argument("-E", "--editSiteAnno", dest = "editSiteAnno", type = "character",
    required = T, help = "editing site annotation file, bed format with additional cols")

parser$add_argument("--annoCols", dest = "annoCols", type = "character",
    required = T, help = "colnames of the editing site annotation file, separted by coma. e.g. chr,start")

parser$add_argument("-N", "--ncRNA", dest = "ncRNA", type = "character",
    required = T, help = "non-conding RNA annotation. 6 col BED format")

parser$add_argument("-O", "--fout", dest = "fout", type = "character",
    required = T, help = "out file")

parser$add_argument("-c", "--core", dest = "core", type = "integer",
    default = 1, help = "Number of core to use, default 1.")

args = parser$parse_args()

coloc.f = args$coloc
anSites = args$editSiteAnno
anSiteCols = args$annoCols
# anSiteCols = c("chr", "start", "end", "edit", "score", "strand", "genefeature", "genename", "genetype", "repeat", "db")
ncRNA.f = args$ncRNA
fout = args$fout
ncore = args$core
cat(paste0("\n\nUsing ", ncore, " cores.\n\n"))



#---------------------------Library---------------------------------



suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(purrr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(glue))

if (ncore > 1) {
    suppressMessages(library(furrr))
    plan(multisession, workers = ncore)
}


#---------------------------Funcs---------------------------------

stringToGr = function(v) {
    #' @description : convert strings into standard format such that GRanges
    #' can convert a character vector into GRange using as(x, 'GRanges')
    #' @param v : character vector
    #' @return  : fixed charactor vector, ready for as(x, 'GRanges')

    pmode1 = function(v) v
    pmode2 = function(v) {
        v = str_split(v, ":")[[1]]
        v = glue("chr{v[[1]]}:{v[[2]]}-{v[[2]]}:*")
        return(v)
    }
    pmode3 = function(v) {
        v = str_split(v, "_")[[1]]
        v = glue("{v[[1]]}:{v[[2]]}-{v[[2]]}:*")
        return(v)
    }
    funcs = list("1" = pmode1, "2" = pmode2, "3" = pmode3)

    pattern = "(^chr[XYM\\d]+:\\d+\\-\\d+:[\\+\\-\\*])|(^\\d+:\\d+:)|(^chr[XYM\\d]+_\\d+_[\\+\\-\\*\\.])"
    pmodes = str_match_all(v, pattern) %>% map_chr(~which(!is.na(.x[2:4])))

    fixed = map2_chr(v, pmodes, ~funcs[[.y]](.x))

    return(fixed)
}

#---------------------------Load---------------------------------

coloc = fread(coloc.f)

coloc[, `:=`(
    PP.H0.abf = NULL, PP.H1.abf = NULL, PP.H2.abf = NULL, PP.H3.abf = NULL
)]


anSiteCols = str_split(anSiteCols, ",")  %>% unlist
anSites = fread(anSites, col.names = anSiteCols)
anSites = anSites[
    , .(pheno1 = paste(chr, end, strand, edit, sep = "_"),
        genefeature, genename, genetype, `repeat`, db)
]
# only keep annotation for editing sites in coloc results
cat(glue("\n\nAll edQTL sites in annotation? {all(coloc$pheno1 %in% anSites$pheno1)}\n\n"))
anSites = anSites[pheno1 %in% coloc$pheno1]


# Noncoding RNA annotation
ncRNA = fread(ncRNA.f, col.names = c("chr", "start","end", "name", "score", "strand"))
ncRNA = makeGRangesFromDataFrame(ncRNA, keep.extra.columns = T)


#---------------------------Main---------------------------------

# right outer join to get annotation, note 1 site may have multiple anno
cat("\n\nAnnotate editing site with gene features.\n")
coloc = coloc[anSites, on = "pheno1"]

### TESTING ONLY >>>
# coloc = head(coloc, 10)

### <<< TESTING ONLY

# determine if top SNP in overlap
cat(glue("\n\nNumber of rows in coloc summary: {nrow(coloc)}.\n\n"))

cat("\n\nCheck if top SNP from each trait is in overlapping region...\n")

pheno1.gr = as(coloc$pheno1  %>% stringToGr, "GRanges")
p1window.gr = as(coloc$p1window, "GRanges")
p2window.gr = as(coloc$p2window, "GRanges")
olapwindow.gr = as(coloc$olapwindow, "GRanges")
p1topsnp.gr = as(coloc$p1topsnp  %>% stringToGr, "GRanges")
p2topsnp.gr = as(coloc$p2topsnp  %>% stringToGr, "GRanges")

p1topin = future_map_chr(seq_along(p1topsnp.gr),
    ~countOverlaps(p1topsnp.gr[.x], olapwindow.gr[.x],
                   minoverlap = 1, type = "within"))

p2topin = future_map_chr(seq_along(p2topsnp.gr),
    ~countOverlaps(p2topsnp.gr[.x], olapwindow.gr[.x],
                   minoverlap = 1, type = "within"))

edsitein = future_map_chr(seq_along(pheno1.gr),
    ~countOverlaps(pheno1.gr[.x], olapwindow.gr[.x],
                   minoverlap = 1, type = "within"))

# pseudo binary code, where,
# `11` = top snp of both phenotypes are in overlap window
# `10` = top snp of trait 1 in overlap window
# `01` = top snp of trait 2 in overlap window
# `00` = neither traits have top snp in overlap window
# overlap window = the overlap between trait 2 qtl window and coloc window
# (which is defined using trait 1 qtl)
topsnpsin = paste(p1topin, p2topin, sep = "")

coloc[, `:=`(
    topsnpsin = topsnpsin,
    edsitein = edsitein
)]

# annotate whether editing site is in NCRNA 1: yes, 0: no
cat("\n\nCheck if editing site is in non-coding RNA.\n")
coloc[, edsiteinNC := as.integer(
    countOverlaps(pheno1.gr, ncRNA, minoverlap = 1, type = "within") > 0
    ) %>% as.character]



#---------------------------Writeout---------------------------------
cat('\n\nWrite out...')
fwrite(coloc, fout, sep = "\t")
cat('Done.\n\n')



