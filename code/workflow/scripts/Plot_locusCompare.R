#---------------------------------------------------------------
#
#' @author      : Chao Dai
#' @Description : Plot locusCompare plots for coloc results
#' @note        : locuscomparer requires internet connection!
#' 
#---------------------------------------------------------------


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
library(locuscomparer)
suppressMessages(library(gridExtra))

#---------------------------------------------------------------
#                 SET UP RUN MODE
#---------------------------------------------------------------

# RUNMODE:
# 1: snakemake
# 2: interactive or jupyter
# Note interactive() returns FALSE when run on jupyter
if (exists("snakemake")) RUNMODE = 1 else RUNMODE = 2

if (RUNMODE == 1) {
    print("### Running in snakemake script mode.")

    INPUT.COLOC.FULL.SUM = snakemake@input[[1]] # collected coloc summary
    OUTPUT = snakemake@output[[1]]

} else {
    INPUT.COLOC.FULL.SUM = "Results/hs38/coloc/RA/coloc_loci_summary.tsv"

}

print(paste0("### INPUT: ", INPUT.COLOC.FULL.SUM))
#study = str_remove_all(INPUT.COLOC.FULL.SUM, "(Results/hs38/coloc/)|(/coloc_loci_summary.tsv)")
study = str_match(INPUT.COLOC.FULL.SUM, "coloc\\/(.+)\\/coloc_loci")[1, 2]
print(paste0("### GWAS study: ", study))
coloc.sum = fread(INPUT.COLOC.FULL.SUM)
coloc.sum[, study := study]

plotLocusCompare = function(dt, STUDY, LOCUS, PID) {
    #' @note : This function requires internet connection.
    
    gwas = dt[study == STUDY & locus == LOCUS & pid == PID, .(rsid = SNP, pval = p.gwas)]
    qtl = dt[study == STUDY & locus == LOCUS & pid == PID, .(rsid = SNP, pval = p.qtl)]
    f.gwas = "gwas.plotLocusCompare.tmp.txt"
    f.qtl = "qtl.plotLocusCompare.tmp.txt"
    write_tsv(gwas, f.gwas)
    write_tsv(qtl, f.qtl)

    g = locuscompare(f.gwas, f.qtl, 
                     title1 = "GWAS", 
                     title2 = "edQTL", 
                     population = "AFR",
                     genome = "hg38")
    return(g)
    
}

# get locus-pid pairs
locus_pid = coloc.sum[, .(locus, pid)] %>% unique
plot_titles = paste(toupper(study), locus_pid$locus, locus_pid$pid, sep=" | ") %>% 
                str_remove('_[\\+\\-]_[AT]_[CG]')

# plot for each locus-pid pair
ggs = map2(locus_pid$locus, locus_pid$pid,
           ~ plotLocusCompare(coloc.sum, study, .x, .y)
           )

ggs.togo = marrangeGrob(grobs = ggs, nrow=1, ncol=1, top = quote(plot_titles[g]))

ggsave(OUTPUT, ggs.togo, width=8, height=5)










