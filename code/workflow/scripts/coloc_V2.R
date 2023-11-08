#---------------------------------------------------------------
#
#' @author      : Chao Dai
#' @Description : This scripted is intended to run by chromosome. It
#'                handles both quant and case-control types of GWAS,
#'                but only quant type of QTL.
#'                Flow chart demonstrating procedures:
#'                https://drive.google.com/file/d/1kAhHic0FqP2OPrl4STaELG2OdcugsC_v/view?usp=sharing
#' @Description : Key inputs GWAS and QTL summary stats must be
#'                processed and formatted in a particular format.
#'                In addition, GWAS config files must be available
#'                to provide key parameters for coloc.
#
#---------------------------------------------------------------

library(tidyverse)
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
library(coloc)
library(yaml)
library(tictoc)

#---------------------------------------------------------------
#                 SET UP RUN MODE
#---------------------------------------------------------------

# RUNMODE:
# 1: snakemake
# 2: interactive or jupyter
# Note interactive() returns FALSE when run on jupyter
if (exists("snakemake")) RUNMODE = 1 else RUNMODE = 2

if (RUNMODE == 1) {
    INPUT.GWAS.SUMSTATS = snakemake@input[['GWAS_sumstats']]
    INPUT.GWAS.SUMSTATS.COLS = snakemake@input[['GWAS_sumstats_cols']]
    INPUT.GWAS.HITS = snakemake@input[['GWAS_hits']]
    INPUT.GWAS.CONFIG = snakemake@input[['GWAS_config']]
    INPUT.QTL.SUMSTATS = snakemake@input[['QTL_sumstats']]
    INPUT.QTL.HITS = snakemake@input[['QTL_hits']]
    OUTPUT.SUMMARY = snakemake@params[['coloc_summary']] # actually in params
    QTL_th = as.numeric(snakemake@params[['QTL_threshold']])
    earlyStop = snakemake@params[['earlystop']]
    Chr = snakemake@wildcards[['chr']]
    GWAS_trait = snakemake@wildcards[['GWAS_trait']]
} else {
    INPUT.GWAS.SUMSTATS = 'Data/GWAS/MS_sumstats_hg38.tsv.gz'
    INPUT.GWAS.SUMSTATS.COLS = 'Data/GWAS/MS_sumstats_hg38.colnames.txt'
    INPUT.GWAS.HITS = 'Data/GWAS/MS_Hits_1e-7_hg38.tsv'
    INPUT.GWAS.CONFIG = "Data/GWAS/gwas_config.yaml"
    INPUT.QTL.SUMSTATS = 'Results/hs38/edQTL_FixedPCs/ranknorm/cis_1000000/chr20/Nominal.txt.gz'
    INPUT.QTL.HITS = 'Results/hs38/edQTL_FixedPCs/ranknorm/cis_10000/permutations_pass_all_chr.AddedQvalue.txt.gz'
    OUTPUT.SUMMARY = 'Results/hs38/TESTCOLOC/test.summary.tsv'
    QTL_th = 0.1
    earlyStop = 'Results/hs38/TESTCOLOC/test.earlystop'
    Chr = "chr20"
    GWAS_trait = "MS"
}

print("### snakemake passed-in strings: ")
walk2(.x = c(rep("input", 6), rep("output", 1), "QTL_FDR_threshold", "chr", "GWAS trait"),
      .y = c(INPUT.GWAS.SUMSTATS, INPUT.GWAS.SUMSTATS.COLS, INPUT.GWAS.HITS,
             INPUT.GWAS.CONFIG, INPUT.QTL.SUMSTATS, INPUT.QTL.HITS,
             OUTPUT.SUMMARY, QTL_th, Chr, GWAS_trait),
      .f = ~ paste(.x, .y, sep = " : ") %>% print
    )

print(paste0("### Work on GWAS trait: ", GWAS_trait, ", on chromosome: ", Chr))


#---------------------------------------------------------------
#                  USEFUL VARIABLES
#---------------------------------------------------------------

# Change column specs if format changed
# this is based on QTLtools v1.3.1 output
QTL_NOMINAL_COLS = c("phenotype_id", "phenotype_chr", "phenotype_start",
                     "phenotype_end", "phenotype_strand", "num_variant_tested",
                     "variant_dist", "genotype_id", "genotype_chr",
                     "genotype_start", "genotype_end","pval_nominal",
                     "pval_r2", "slope", "top_variant_flag")
QTL_NOMINAL_COLCLASSES = c("character", "character", "integer",
                           "integer", "character", "integer",
                           "integer", "character", "character",
                           "integer", "integer", "double",
                           "double", "double", "integer")
cis_window = 10000 # 10kb cis window for QTL hits
GWAS_config = read_yaml(INPUT.GWAS.CONFIG)
GWAS_config = GWAS_config$GWAS[[GWAS_trait]]


#---------------------------------------------------------------
#                  USEFUL FUNCTIONS
#---------------------------------------------------------------


ReadQtlNominalByRegion = function(x, filename, COLNAMES, COLCLASSES) {
  #' @method           : Use tabix and fread to load a specific region of summary
  #'                     statistics of QTL nominal pass.
  #' @param x          : a single genomic range that specifies regions to read from.
  #' @param filename   : file name of the qtl nominal pass file.
  #' @param COLNAMES   : a vector specifying colnames.
  #' @param COLCLASSES : a vector specifying column datatypes.
  #' @return           : data.table for QTL nominal SNPs subset by x

  region = as.character(x)
  stats = fread(
    cmd = paste0("tabix ", filename, " ", region),
    header = F,
    col.names = COLNAMES,
    colClasses = COLCLASSES
  ) %>%
    unique

  return(stats)
}


RunColoc = function(gwas, qtl,
                    gwas_type, qtl_type = "quant",
                    gwas_N, qtl_N = 85,
                    gwas_s = NA, gwas_sdY = NA,
                    qtl_sdY = 1) {
  #' @description     : Run coloc on 1 coloc window with 1 set of GWAS SNPs
  #'                    and 1 set of phenotype QTL snps.
  #' @note            : When type is "quant", coloc requires sdY, if sdY is not 1.
  #'                    Or estimate sdY using MAF and N. When type is "cc", coloc
  #'                    requires the "s" parameter.
  #' @param gwas      : GRange object of GWAS SNPs with sumstats.
  #' @param qtl       : GRange object of QTL SNPs with sumstats of 1 phenotype.
  #' @param gwas_type : Must be either "quant" or "cc".
  #' @param qtl_type  : Must be either "quant" or "cc", default to "cc".
  #' @param gwas_N    : N, total sample numbers in GWAS.
  #' @param qtl_N     : N, total sample numbers in QTL.
  #' @param gwas_s    : The s parameter required by coloc for case-control(cc) type.
  #' @param gwas_sdY  : sdY parameter for GWAS.
  #' @param qtl_sdY   : sdY parameter for QTL. set to 1, since its normalized.
  #' @return          : A two level nested list of coloc.abf result.


  gwas_mcolnames = mcols(gwas) %>% names
  qtl_mcolnames = mcols(qtl) %>% names

  # first convert GRange to data.table
  gwas = as.data.table(gwas) %>% .[, -c("width", "strand")]
  qtl = as.data.table(qtl) %>% .[, -c("width", "strand")]

  # some SNPs are actually INDELs, thus 1 snp_id may have multiple genotype_id
  # only keep 1 row per snp_id (use the one with the lowest pval_nominal)
  gwas[, prank := rank(pval, ties.method = "first"), by = snp_id]
  gwas = gwas[prank == 1 & se != 0] # snps with se==0 will throw error
  qtl[, prank := rank(pval_nominal, ties.method = "first"), by = snp_id]
  qtl = qtl[prank == 1]

  # compute z and varbeta for QTL sumstats
  qtl[, z := qnorm(pval_nominal / 2, lower.tail = F)] # compute z
  qtl[, varbeta := (slope / z)^2] # compute variance
  qtl = qtl[!is.na(varbeta)] # remove Nan

  if (gwas_type == "cc") {
    # Case-control mode -------------------------------------------------------

    if (is.na(gwas_s)) stop("Error. Case-control GWAS needs the S paramter specified.")

    ds1 = list(
      beta = gwas$beta,
      varbeta = (gwas$se)^2,
      type = gwas_type,
      s = gwas_s,
      snp = gwas$snp_id,
      position = gwas$start,
      N = gwas_N)

    ds2 = list(
      beta = qtl$slope,
      varbeta = qtl$varbeta,
      type = qtl_type,
      sdY = qtl_sdY,
      snp = qtl$snp_id,
      position = qtl$start,
      N = qtl_N)

    res = suppressMessages(coloc.abf(ds1,ds2))

    # case-control ends next }
  } else if (gwas_type == "quant") {
    # Quant mode --------------------------------------------------------------
    if (!is.na(gwas_sdY)) { # has sdY value specified
      ds1 = list(
        beta = gwas$beta,
        varbeta = (gwas$se)^2,
        type = gwas_type,
        sdY = gwas_sdY,
        snp = gwas$snp_id,
        position = gwas$start,
        N = gwas_N)
    } else if ("maf" %in% names(gwas)) { # has maf column
      ds1 = list(
        beta = gwas$beta,
        varbeta = (gwas$se)^2,
        type = gwas_type,
        MAF = gwas$maf,
        snp = gwas$snp_id,
        position = gwas$start,
        N = gwas_N)
    } else {
      stop("Error. GWAS needs either the 'maf' column or has 'sdY' specified.")
    }

    ds2 = list(
      beta = qtl$slope,
      varbeta = qtl$varbeta,
      type = qtl_type,
      sdY = qtl_sdY,
      snp = qtl$snp_id,
      position = qtl$start,
      N = qtl_N)

    res = suppressMessages(coloc.abf(ds1,ds2))

  } # quant mode ends.

  return(res)
} # function ends.


hush = function(code) {
    # shut up coloc. credit: stackoverflow
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
}


#---------------------------------------------------------------
#              STEP1. GET COLOC FOCUS WINDOWS
#---------------------------------------------------------------
#' @description Coloc focus windows are GWAS hit loci that contains
#'              at least 1 edQTL hit.



# Get GWAS hit loci -------------------------------------------------------

print(paste0("### Working on chromosome ", Chr))
print("### Get GWAS hit loci")

gwas_hits = fread(INPUT.GWAS.HITS) %>%
    .[seqnames %in% Chr] %>%
    unique

# in case of 0 hits in gwas
if (nrow(gwas_hits) == 0){
    readr::write_lines(paste0(Chr, " has 0 GWAS hits. Quit script."), earlyStop)
    quit(save = "no")
}

# Convert to GRange object
gwas_hits = makeGRangesFromDataFrame(gwas_hits, keep.extra.columns = T)



# Remove HLA region -------------------------------------------------------

if (Chr == "chr6") {
    print(paste0("Remove region if overlapping at least 5kb with HLA regions: chr6:24999772-35032223..."))
    hla.gr = GRanges(seqnames = "chr6", ranges = IRanges(start=24999772, end=35032223))
    hits.gwas = suppressWarnings(
      subsetByOverlaps(gwas_hits,
                       ranges = hla.gr,
                       minoverlap = 5000,
                       invert = T)
      )
}

# in case of 0 hits in gwas
if (length(gwas_hits) == 0){
    readr::write_lines(paste0(Chr, " has 0 GWAS hits. Quit script."), earlyStop)
    quit(save = "no")
}

print(paste0("Found ", length(gwas_hits), " GWAS hit loci on chromosome ", Chr, "."))





# Get QTL hits ------------------------------------------------------------

print(paste0("### Get edQTL permutation pass"))
qtl_hits = read_delim(INPUT.QTL.HITS,
                      delim = " ",
                      col_types = "cciiciicciiiddddddddd") %>%
  dplyr::filter(q < QTL_th & phenotype_chr %in% !!Chr)

# in case of 0 hits in QTL
if (nrow(qtl_hits) == 0){
  readr::write_lines(paste0(Chr, " has 0 edQTL hits. Quit script."), earlyStop)
  quit(save = "no")
}

print(paste0("Loaded ",
             nrow(qtl_hits),
             " edQTL permutation pass hits with FDR<",
             QTL_th," on ", Chr, "."))

qtl_hits_gr = makeGRangesFromDataFrame(qtl_hits,
                                    keep.extra.columns = T,
                                    seqnames.field = "phenotype_chr",
                                    start.field = "phenotype_start",
                                    end.field = "phenotype_start")




# Determine coloc focus windows -------------------------------------------
#' @description GWAS hit locus must have >=1 edQTL hits

coloc_focus_gr = subsetByOverlaps(gwas_hits, qtl_hits_gr, minoverlap = 1) %>% unique

# if no coloc focus window exist, quit.
if (length(coloc_focus_gr) == 0) {
    readr::write_lines(
      paste0(Chr,
             " has 0 overlapping hits between GWAS and edQTL. Quit script."),
      earlyStop)
    quit(save = "no")
}

coloc_focus_ls = split(coloc_focus_gr, 1:length(coloc_focus_gr)) %>% as.list
names(coloc_focus_ls) = as.character(coloc_focus_gr)

print(paste0("Found ", length(coloc_focus_ls), " GWAS hit loci containing edQTL passing FDR."))



#---------------------------------------------------------------
#               STEP2. SUBSET edQTL HITS
#---------------------------------------------------------------
#' @description Only include edQTL hits that overlap with coloc
#'              focus windows. Keep the phenotype_id as key and
#'              later use as filtering criteria in edQTL sumstats.


# resize edQTL to 10kb cis-window
qtl_hits_cis_region = promoters(qtl_hits_gr,
                                upstream = as.integer(cis_window/2),
                                downstream = as.integer(cis_window/2)
                                )

# require at least 10% of cis-window overlaps with coloc focus window
qtl_hits_cis_region = map(coloc_focus_ls,
                          ~subsetByOverlaps(x = qtl_hits_cis_region,
                                            ranges = .x,
                                            minoverlap = as.integer(.1*cis_window)) %>%
                            unique
                          )

# make sure each focus windows have edQTL left
qtl_hits_cis_region = qtl_hits_cis_region[lengths(qtl_hits_cis_region) > 0]

# if no coloc focus window exist, quit.
if (length(qtl_hits_cis_region) == 0) {
  readr::write_lines(paste0(Chr, " has 0 overlapping hits between GWAS and edQTL. Quit script."), earlyStop)
  quit(save = "no")
}

QTL_IDS_IN_FOCUS = purrr::reduce(qtl_hits_cis_region, c) %>%
  .$phenotype_id %>%
  unique
print(paste0("### Found ", length(QTL_IDS_IN_FOCUS), " edQTL qualified for coloc analysis."))



#---------------------------------------------------------------
#         STEP3. GWAS SUMSTATS SNPS FOR COLOC
#---------------------------------------------------------------
#' @description For each coloc focus window, get the set of GWAS
#'              SNPs and summary stats.


# For each focus window, get overlapping GWAS SNPs ------------------------

print(paste0("### Load GWAS summary stats for chromosome ", Chr, "."))
stats.gwas.cols = read_lines(INPUT.GWAS.SUMSTATS.COLS)
stats_gwas = fread(cmd = paste0("tabix ", INPUT.GWAS.SUMSTATS, " ", Chr),
                   col.names = stats.gwas.cols
                   ) %>% unique
print(paste0("Found ", nrow(stats_gwas), " GWAS SNPs on ", Chr, "."))

stats_gwas = makeGRangesFromDataFrame(stats_gwas,
                                      keep.extra.columns = T,
                                      seqnames.field = "chr",
                                      start.field = "pos",
                                      end.field = "pos")

# subset gwas stats to only include SNPs that are in coloc windows
stats_gwas = map(coloc_focus_ls,
                 ~ subsetByOverlaps(stats_gwas, .x, minoverlap = 1) %>%
                   unique
                 )

iwalk(stats_gwas,
      ~print(paste0("Found ", length(.x), " overlapping GWAS SNPs in focus: ", .y))
      )



#---------------------------------------------------------------
#           STEP4. QTL SUMSTATS SNPS FOR COLOC
#---------------------------------------------------------------
#' @description Find a list of list, such that level 1 stores the
#'              coloc focus windows; level 2 stores the corresponding
#'              edQTL's SNPs (sumstats) for each coloc window


# EDITING QTL summary stats (nominal pass results)
print(paste0("### Load edQTL nominal pass SNPs for chromosome ", Chr, "."))


# QTL sumstats for phenotypes that are in coloc windows
print("Include only edQTL within coloc window.")
stats_qtl = map(coloc_focus_ls,
                ~ReadQtlNominalByRegion(x = .x,
                                        filename = INPUT.QTL.SUMSTATS,
                                        COLNAMES = QTL_NOMINAL_COLS,
                                        COLCLASSES = QTL_NOMINAL_COLCLASSES) %>%
                  .[phenotype_id %in% QTL_IDS_IN_FOCUS]
                )

# for each coloc window and each phenotype, subset include only SNPs in coloc window
stats_qtl = map(stats_qtl,
                ~makeGRangesFromDataFrame(.x,
                                          keep.extra.columns = T,
                                          seqnames.field = "genotype_chr",
                                          start.field = "genotype_start",
                                          end.field = "genotype_end"
                                          )
                )
if (all(names(stats_qtl) == names(coloc_focus_ls))) {
  stats_qtl = map2(stats_qtl,
                   coloc_focus_ls,
                   ~subsetByOverlaps(.x, .y, minoverlap = 1)
                   )
}

# split each coloc window into a list by phenotype_id
stats_qtl = map(stats_qtl,
                ~split(.x, .x$phenotype_id) %>% as.list
                )
# create a snp_id column to make ensure stat_qtl and stat_gwas can be joined

print("Include only edQTL nominal SNPs in coloc window. Remove any edQTL with less than 20 shared SNPs.")
stats_qtl = map(stats_qtl,
                ~ map(.x,
                     function(gr) {
                       gr$snp_id = paste(seqnames(gr), start(gr), sep=":")
                       return(gr)
                     }
                     )
                )

# snp_id in each QTL must be shared with stat_gwas SNPs, min shared snps >=20
stats_qtl = map2(.x = stats_qtl,
                  .y = stats_gwas,
                  .f = ~ map(.x,
                             function(gr) {
                               rowidx = gr$snp_id %in% .y$snp_id # QTL's snp_id must be in gwas's snp_id
                               return(gr[rowidx]) # some QTL may have 0 shared snps
                             }
                  )
)
stats_qtl = map(stats_qtl, ~.x[lengths(.x) >= 20])


print("### summary stastatics ready for colocalization analysis: ")
num.qtl.per.locus = map_int(stats_qtl, length)
num.qtlsnp.per.locus = map(stats_qtl,
                           ~purrr::reduce(.x, c) %>%
                             .$snp_id %>%
                             unique %>%
                             length
                           )
num.gwassnp.per.locus = lengths(stats_gwas)
pwalk(.l = list(coloc_focus_ls,
                num.gwassnp.per.locus,
                num.qtl.per.locus,
                num.qtlsnp.per.locus),
      .f = ~ print(paste0(
        "Coloc locus: ", ..1, " ; ",
        "GWAS SNPS: ", ..2, " ; ",
        "Number of edQTL: ", ..3, " ; ",
        "Number of SNPs shared with GWAS across edQTL: ", ..4, "."))
      )




#---------------------------------------------------------------
#                 STEP5. RUN COLOC BY WINDOW BY PHENOTYPE
#---------------------------------------------------------------


# Prepare parameters for coloc --------------------------------------------

gType = GWAS_config$type

print("### Run colocalization: ")
if (gType == "cc") {
  gCaseN = as.integer(GWAS_config$N_case)
  gControlN = as.integer(GWAS_config$N_control)
  gN = gCaseN + gControlN
  gS = gCaseN / gControlN
  gsdY = NA

  print(paste0("GWAS type: ", gType, " , ",
               "Case: ", gCaseN, " , Control: ", gControlN, " , ",
               "Proportion case: ", round(gS,2), "."
               ))

} else if (gType == "quant") {
  gN = GWAS_config$N
  gS = NA
  if (length(GWAS_config$sdY) > 0) {
    gsdY = as.numeric(GWAS_config$sdY)
    print(paste0("GWAS type: ", gType, ", provided sdY: ", gsdY, "."))
  } else if (length(GWAS_config$maf) > 0) {
    gsdY = NA
    gN = GWAS_config$N
    print(paste0("GWAS type: ", gType, ", no sdY, but provided maf."))
    print("coloc will use varbeta and maf to estimate sdY.")
  }
}


# Run coloc.abf -----------------------------------------------------------

results = map2(.x = stats_qtl, # nested list of qtl colocs, level2 is per phenotype
               .y = stats_gwas, # list of gwas snps in locus
               .f = ~ map(.x = .x, # list of SNPs of a phenotype of a locus
                          function(xx) { # xx: qtl nominal SNPs of 1 phenotype
                            res = hush(RunColoc(
                              gwas = .y,
                              qtl = xx,
                              gwas_type = gType,
                              gwas_N = gN,
                              gwas_s = gS,
                              gwas_sdY =gsdY
                            ))}
                          )
               )



# Collect and write out results -------------------------------------------


results.summary = imap_dfr(
  results,
  ~ imap_dfr(.x,
             function(x,y) x$summary %>%
               enframe %>%
               pivot_wider(names_from = name) %>%
               add_column(pid = y, locus = .y)
             )
  )

results.summary = as.data.table(results.summary) %>%
  .[order(-PP.H4.abf,pid),
    .(locus, pid, PP.H4.abf, nsnps)]


print("### Write output to files")
write_tsv(results.summary, OUTPUT.SUMMARY)

print("### Finished ###")

