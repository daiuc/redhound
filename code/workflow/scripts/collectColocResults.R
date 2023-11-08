suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))

# RUNMODE:
# 1: snakemake
# 2: interactive or jupyter
# Note interactive() returns FALSE when run on jupyter
if (exists("snakemake")) RUNMODE = 1 else RUNMODE = 2

if (RUNMODE == 1) {
    INPUT.GWAS.SUMSTATS = snakemake@input[['GWAS_sumstats']]
    INPUT.GWAS.SUMSTATS.COLS = snakemake@input[['GWAS_sumstats_cols']]
    INPUT.dbsnp150 = snakemake@input[['dbsnp150_lookup']]
    OUTPUT.dt = snakemake@output[[1]]
    QTL.SUMSTATS.prefix = snakemake@params[['QTL_sumstats_prefix']]
    INPUT.COLOC.SUMS = list.files(paste(snakemake@params[['COLOC_summary_prefix']], 1:22,sep = "" ),
                                  "summary.tsv",
                                  full.names = T)
    GWAS_trait = snakemake@wildcards[['GWAS_trait']]

} else {
    INPUT.GWAS.SUMSTATS = 'Data/GWAS/Height_sumstats_hg38.tsv.gz'
    INPUT.GWAS.SUMSTATS.COLS = 'Data/GWAS/Height_sumstats_hg38.colnames.txt'
    INPUT.dbsnp150 = "/project2/yangili1/cdai/SNP/snp150_hg38_3cols.txt.gz"
    OUTPUT.dt = "colocResult.test"
    QTL.SUMSTATS.prefix = "Results/hs38/edQTL_FixedPCs/ranknorm/cis_1000000/"
    INPUT.COLOC.SUMS = list.files(paste("Results/hs38/coloc/Height/chr", 1:22,sep = "" ),
                                    "summary.tsv",
                                    full.names = T)
    GWAS_trait = "Height"
}


print("### snakemake passed-in strings: ")
print(paste0("### Collect coloc summary files for ", GWAS_trait, " and edQTL."))
autosome = paste("chr", 1:22, sep="")
print("### Load colocalization summary files: ")
print(INPUT.COLOC.SUMS)

# Useful variables --------------------------------------------------------


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


# Useful functions --------------------------------------------------------


ReadQtlNominalByRegion = function(x, filename, COLNAMES=NA, COLCLASSES=NA) {
    #' @method           : Use tabix and fread to load a specific region of summary
    #'                     statistics of QTL nominal pass.
    #' @param x          : a single genomic range that specifies regions to read from.
    #' @param filename   : file name of the qtl nominal pass file.
    #' @param COLNAMES   : a vector specifying colnames.
    #' @param COLCLASSES : a vector specifying column datatypes.
    #' @return           : data.table for QTL nominal SNPs subset by x

    region = as.character(x)
    print(paste0("tabix ", filename, " ", region))
    stats = fread(
        cmd = paste0("tabix ", filename, " ", region),
        header = F,
        col.names = COLNAMES,
        colClasses = COLCLASSES
    ) %>%
        unique

    return(stats)
}


# Get coloc regions -------------------------------------------------------


coloc_summary = map(INPUT.COLOC.SUMS, fread) %>% rbindlist
coloc_summary[, `:=`(
    chr = factor(str_extract(locus, "chr[0-9]{1,2}"), autosome),
    start = as.integer(str_split(locus, "[:-]", simplify = T)[,2]),
    end = as.integer(str_split(locus, "[:-]", simplify = T)[,3])
    )]

# select loci with PP.H4.abv > 0.8
coloc_summary = coloc_summary[PP.H4.abf > .8][order(PP.H4.abf)]
coloc_summary = coloc_summary[!is.na(start) & !is.na(end)]
# if no coloc focus window exist, quit.
if (nrow(coloc_summary) == 0) {
    print("No colocalizing locus found at PP.H4.abf > .8")
    quit(save = "no")
}


print(paste0("### Found ", nrow(coloc_summary), " colocalized loci + pid passing PP.H4.abv > 0.8"))
coloc_regions = makeGRangesFromDataFrame(coloc_summary[, .(chr, start, end)],
                                         keep.extra.columns = F) %>% unique

print(paste("### Colocalized loci:", sort(coloc_regions)))
coloc_regions = split(coloc_regions, 1:length(coloc_regions)) %>% as.list
names(coloc_regions) = map_chr(coloc_regions, as.character)




# Read in summary stats ---------------------------------------------------

print(paste0("Load GWAS sumstats for the regions."))
GWAS.SUMSTATS.colnames = read_lines(INPUT.GWAS.SUMSTATS.COLS)
stats_gwas = map(coloc_regions,
                 ~ ReadQtlNominalByRegion(.x,
                                          INPUT.GWAS.SUMSTATS,
                                          GWAS.SUMSTATS.colnames
                                          )
                 )

print(paste0("Load QTL sumstats for the regions."))
qtl.sum.files = map(coloc_regions,
                ~ paste(QTL.SUMSTATS.prefix, seqnames(.x), '/Nominal.txt.gz', sep="")
                    )
QTL_IDS_IN_COLOC = coloc_summary$pid %>% unique
stats_qtl = map2(coloc_regions, qtl.sum.files,
                ~ ReadQtlNominalByRegion(.x,
                                         .y,
                                         QTL_NOMINAL_COLS,
                                         QTL_NOMINAL_COLCLASSES
                                         ) %>%
                    .[phenotype_id %in% QTL_IDS_IN_COLOC,
                      .(phenotype_id, pval_nominal, slope, top_variant_flag,
                        genotype_chr, genotype_start, genotype_id,
                        snp_id = paste(genotype_chr, genotype_start, sep=":"))
                      ]
                )

print("Create flat datframes for gwas and qtl for the snps in coloc regions.")
shared_snps = map2(stats_gwas, stats_qtl, ~intersect(.x$snp_id, .y$snp_id))

# keep only shared SNPs
stats_gwas = map2(stats_gwas, shared_snps, ~.x[snp_id %in% .y])
stats_qtl = map2(stats_qtl, shared_snps, ~.x[snp_id %in% .y])

# flatten nested lists
stats_gwas = imap(stats_gwas, ~.x[,locus := .y]) %>% rbindlist
stats_qtl = imap(stats_qtl, ~.x[, locus := .y]) %>% rbindlist



# Annotate with rs id -----------------------------------------------------
print("### Annotate SNPs with 'rs' SNP ids using dbSNP150.")

dbsnp_regions = names(coloc_regions)
dbsnp_regions = paste("tabix", INPUT.dbsnp150, dbsnp_regions, sep=" ")
dbsnp150 = map(dbsnp_regions, ~fread(cmd = .x, header = F, sep = "\t")) %>% rbindlist
names(dbsnp150) = c("CHR","BP", "SNP")

stats_gwas = dbsnp150[stats_gwas, on=c(CHR = "chr", BP = "pos")]
stats_gwas[, SNP := if_else(is.na(SNP), snp_id, SNP)]

stats_qtl = dbsnp150[stats_qtl, on=c(CHR = "genotype_chr", BP = "genotype_start")]
stats_qtl[, SNP := if_else(is.na(SNP), snp_id, SNP)]


# join GWAS snps and coloc locus to bring in PP.H4, and p.gwas
a = left_join(stats_gwas[, .(locus, CHR, BP, SNP, p = pval)],
          coloc_summary,
          by = "locus")

# then join with QTL to bring together p.qtl
dt = inner_join(a, stats_qtl[, .(locus, phenotype_id, SNP, p = pval_nominal)], suffix = c(".gwas", ".qtl"),
           by=c("locus", pid = "phenotype_id", "SNP"),)
dt[, `:=`(
    CHR = str_remove(CHR, "chr") %>% as.integer
)]

write_tsv(dt, OUTPUT.dt)

# plot --------------------------------------------------------------------
#
# source("workflow/Scripts/qqman/R/ggColocManhattan.R")
#
#
# ggColocManhattan(summ.1 = dt[locus == "chr1:155487637-156280569",
#                              .(CHR, BP, SNP, P=p.gwas)],
#                  summ.2 = dt[locus == "chr1:155487637-156280569" & pid == "chr1_155444591_+_A_G",
#                              .(CHR, BP, SNP, P=p.qtl)],
#                  PP4 = .803,
#                  ) +
#     scale_x_continuous(labels = number_format(scale = 1e-6))



