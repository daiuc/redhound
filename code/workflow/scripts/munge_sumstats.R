## Only run this script interactively!!!
## Recommend running each GWAS with a fresh session to avoid memory overflow

library(tidyverse)
library(data.table)
library(furrr)
library(tictoc)


NCore = as.integer(system2("nproc", stdout = T))
NCore = min(NCore, 20)
print(paste0("Number of Cores: ", NCore))

# chain = import.chain("/home/chaodai/software/liftover/hg19ToHg38.over.chain")
GWAS.dir.prefix = "/project2/yangili1/cdai/aicher/code/chRNA/Data/GWAS/gwas_download/"
GWAS.out.prefix = "/project2/yangili1/cdai/aicher/code/chRNA/Data/GWAS/"
autosomes = paste("chr", 1:22, sep="")



# -------------------------------------------------------------------------









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                Manually Munge Summary Stats!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#-----------------------------------------------------------
#               Height Height Height
#-----------------------------------------------------------

Height = fread(paste0(GWAS.dir.prefix, "Locke_height_UKBiobank_2018.txt.gz"),
               colClasses = c(rep("integer", 2), rep("character", 3,), rep("double", 5)))

Height[, `:=`(
    chr = factor(paste("chr", CHR, sep = ""), levels = autosomes),
    BEDstart = as.integer(POS) - 1,
    BEDend = as.integer(POS),
    snp_id = paste("chr", CHR, ":", POS, sep=""),
    beta = BETA,
    se = SE,
    pval = P,
    n = N,
    maf = Freq_Tested_Allele_in_HRS,
    snp_ori = SNP# original
)]

# select columns, make sure to sort rows
Height_cols = c("chr", "BEDstart", "BEDend", "snp_id", "beta", "se", "pval", "n", "maf", "snp_ori")
Height = Height[chr %in% autosomes, ..Height_cols][order(chr, BEDstart, BEDend)]


# write 3 column bed file for coordinates
write_tsv(Height[, .(chr, BEDstart, BEDend, snp_id)],
          paste0(GWAS.out.prefix, "Height_stats_coordinates_hg19.bed"),
          col_names = F)


# Liftover from hg19 to hg38 ----------------------------------------------

#deleted.snps = fread(cmd = "awk '$1 !~ /#Dele/' /project2/yangili1/cdai/aicher/code/chRNA/Data/GWAS/Height.liftover.unmapped")

# do it in shell in the GWAS dir
# ./hg19Tohg38.sh my_hg19.bed my_hg38.bed my.unmapped


# Update hg19 cordinates in data.table to hg38 ----------------------------
Height.hg38.bed = fread(paste0(GWAS.out.prefix, "Height_stats_coordinates_hg38.bed"),
                        col.names = c("chr", "BEDstart", "BEDend", "snp_id"),
                        colClasses = c("character", rep("integer", 2), "character")) %>% unique
Height.hg38.bed = Height.hg38.bed[chr %in% autosomes] # in case some hg19 coords mapped to scaffolds
Height.hg38.bed[, chr := factor(chr, autosomes)]


# do inner join to bring in summary stats
Height38 = Height.hg38.bed[Height, on = "snp_id", nomatch = NULL][!is.na(chr)]
# remove hg19 coordinates, and old snp_id based on hg19 coords
Height38[, `:=`(i.chr = NULL, i.BEDstart = NULL, i.BEDend = NULL, snp_id = NULL)]
# reconstruct snp_id using hg38 chr and end position.
Height38[, snp_id := paste(chr, BEDend, sep=":")]
# column order, and ensure row sorted
Height38 = Height38[, ..Height_cols][order(chr, BEDstart, BEDend),] %>% unique

# now saving it to non-bed format. pos is using 1-indexed
Height38 = Height38[, .(chr, pos = BEDend, snp_id, beta, se, pval, n, maf, snp_ori)]

# Finally! Write munged summary stats -------------------------------------

# write munged stats, column names in separate file
write_lines(names(Height38), paste0(GWAS.out.prefix, "Height_sumstats_hg38.colnames.txt"))
write_tsv(Height38,
          file = paste0(GWAS.out.prefix, "Height_sumstats_hg38.tsv"),
          col_names = F)



# Don't forget to bgzip and tabix in shell --------------------------------
# Also safe to delete the intermediary files: *_hg19.bed, *_hg38.bed, but keep the *.liftover.unmapped file.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#-----------------------------------------------------------
#                       IBD IBD IBD
#-----------------------------------------------------------

IBD = fread(paste0(GWAS.dir.prefix, "ibd_build37_59957_20161107.txt.gz"),
               colClasses = c(rep("character", 3), rep("double", 3,),
                              "character", rep("double", 8)))

# split Markername to dt of chr and pos
chr_pos = IBD[, .(MarkerName)
              ][,str_split(MarkerName, "[:_]", simplify = T) %>%
                    as.data.frame %>% as.list %>% .[1:2]
                ]
chr_pos = chr_pos[, `:=`(chr = paste("chr", V1, sep=""), pos = V2)][, .(chr, pos)]


IBD[, `:=`(
    beta = Effect,
    se = StdErr,
    pval = `P.value`,
    snp_ori = MarkerName # original
)]

IBD = cbind(chr_pos, IBD)
IBD[, `:=`(
    chr = factor(chr, autosomes),
    snp_id = paste(chr, pos, sep=":"),
    BEDstart = as.integer(pos) - 1,
    BEDend = as.integer(pos)
)]

# select columns, make sure to sort rows
IBD_cols = c("chr", "BEDstart", "BEDend", "snp_id", "beta", "se", "pval", "snp_ori")
IBD = IBD[chr %in% autosomes, ..IBD_cols][order(chr, BEDstart, BEDend)]


# write 3 column bed file for coordinates
write_tsv(IBD[, .(chr, BEDstart, BEDend, snp_id)],
          paste0(GWAS.out.prefix, "IBD_stats_coordinates_hg19.bed"),
          col_names = F)


# Liftover from hg19 to hg38 ----------------------------------------------

#deleted.snps = fread(cmd = "awk '$1 !~ /#Dele/' /project2/yangili1/cdai/aicher/code/chRNA/Data/GWAS/IBD.liftover.unmapped")

# do it in shell in the GWAS dir
# ./hg19Tohg38.sh my_hg19.bed my_hg38.bed my.unmapped


# Update hg19 cordinates in data.table to hg38 ----------------------------
IBD.hg38.bed = fread(paste0(GWAS.out.prefix, "IBD_stats_coordinates_hg38.bed"),
                        col.names = c("chr", "BEDstart", "BEDend", "snp_id"),
                        colClasses = c("character", rep("integer", 2), "character")) %>% unique
IBD.hg38.bed = IBD.hg38.bed[chr %in% autosomes] # in case some hg19 coords may be lifted to scaffolds
IBD.hg38.bed[, chr := factor(chr, autosomes)]

# do inner join to bring in summary stats (data.table does right outer join)
IBD38 = IBD.hg38.bed[IBD, on = "snp_id", nomatch = NULL][!is.na(chr)]
# remove hg19 coordinates, and old snp_id based on hg19 coords
IBD38[, `:=`(i.chr = NULL, i.BEDstart = NULL, i.BEDend = NULL, snp_id = NULL)]
# reconstruct snp_id using hg38 chr and end position.
IBD38[, snp_id := paste(chr, BEDend, sep=":")]
# column order, and ensure row sorted
IBD38 = IBD38[, ..IBD_cols][order(chr, BEDstart, BEDend),] %>% unique

# now saving it to non-bed format. pos is using 1-indexed
IBD38 = IBD38[, .(chr, pos = BEDend, snp_id, beta, se, pval, snp_ori)]

# Finally! Write munged summary stats -------------------------------------

# write munged stats, column names in separate file
write_lines(names(IBD38), paste0(GWAS.out.prefix, "IBD_sumstats_hg38.colnames.txt"))
write_tsv(IBD38,
          file = paste0(GWAS.out.prefix, "IBD_sumstats_hg38.tsv"),
          col_names = F)



# Don't forget to bgzip and tabix in shell --------------------------------
# Also safe to delete the intermediary files: *_hg19.bed, *_hg38.bed, but keep the *.liftover.unmapped file.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#-----------------------------------------------------------
#                       RA RA RA
#-----------------------------------------------------------

# RA = fread(paste0(GWAS.dir.prefix, "Trans_all_auto-10-2021.txt.gz"),
#            colClasses = c(rep("character", 4)))

RA = read_table(paste0(GWAS.dir.prefix, "Trans_all_auto-10-2021.txt.gz"),
                col_types = "cddd") %>% as.data.table

# split Markername to dt of chr and pos
chr_pos = RA[, .(SNP)
][,str_split(SNP, "[:_]", simplify = T) %>%
      as.data.frame %>% as.list %>% .[1:2]
]
chr_pos = chr_pos[, `:=`(chr = paste("chr", V1, sep=""), pos = V2)][, .(chr, pos)]


RA[, `:=`(
    beta = Beta,
    se = SE,
    pval = Pval,
    snp_ori = SNP # original
)]

RA = cbind(chr_pos, RA)
RA[, `:=`(
    chr = factor(chr, autosomes),
    snp_id = paste(chr, pos, sep=":"),
    BEDstart = as.integer(pos) - 1,
    BEDend = as.integer(pos)
)]

# select columns, make sure to sort rows
RA_cols = c("chr", "BEDstart", "BEDend", "snp_id", "beta", "se", "pval", "snp_ori")
RA = RA[chr %in% autosomes, ..RA_cols][order(chr, BEDstart, BEDend)]


# write 3 column bed file for coordinates
write_tsv(RA[, .(chr, BEDstart, BEDend, snp_id)],
          paste0(GWAS.out.prefix, "RA_stats_coordinates_hg19.bed"),
          col_names = F)


# Liftover from hg19 to hg38 ----------------------------------------------

#deleted.snps = fread(cmd = "awk '$1 !~ /#Dele/' /project2/yangili1/cdai/aicher/code/chRNA/Data/GWAS/RA.liftover.unmapped")

# do it in shell in the GWAS dir
# ./hg19Tohg38.sh my_hg19.bed my_hg38.bed my.unmapped


# Update hg19 cordinates in data.table to hg38 ----------------------------
RA.hg38.bed = fread(paste0(GWAS.out.prefix, "RA_stats_coordinates_hg38.bed"),
                     col.names = c("chr", "BEDstart", "BEDend", "snp_id"),
                     colClasses = c("character", rep("integer", 2), "character")) %>% unique
RA.hg38.bed = RA.hg38.bed[chr %in% autosomes] # in case some hg19 coords may be lifted to scaffolds
RA.hg38.bed[, chr := factor(chr, autosomes)]

# do inner join to bring in summary stats (data.table does right outer join)
RA38 = RA.hg38.bed[RA, on = "snp_id", nomatch = NULL][!is.na(chr)]
# remove hg19 coordinates, and old snp_id based on hg19 coords
RA38[, `:=`(i.chr = NULL, i.BEDstart = NULL, i.BEDend = NULL, snp_id = NULL)]
# reconstruct snp_id using hg38 chr and end position.
RA38[, snp_id := paste(chr, BEDend, sep=":")]
# column order, and ensure row sorted
RA38 = RA38[, ..RA_cols][order(chr, BEDstart, BEDend),] %>% unique

# now saving it to non-bed format. pos is using 1-indexed
RA38 = RA38[, .(chr, pos = BEDend, snp_id, beta, se, pval, snp_ori)]

# Finally! Write munged summary stats -------------------------------------

# write munged stats, column names in separate file
write_lines(names(RA38), paste0(GWAS.out.prefix, "RA_sumstats_hg38.colnames.txt"))
write_tsv(RA38,
          file = paste0(GWAS.out.prefix, "RA_sumstats_hg38.tsv"),
          col_names = F)



# Don't forget to bgzip and tabix in shell --------------------------------
# Also safe to delete the intermediary files: *_hg19.bed, *_hg38.bed, but keep the *.liftover.unmapped file.


#-----------------------------------------------------------
#                       MS MS MS
#-----------------------------------------------------------

MS = fread(paste0(GWAS.dir.prefix, "MS_discovery_metav3.0.meta.gz"),
           colClasses = c(rep("integer", 2), rep("character", 3),
                          "integer", rep("double",2)))

# check if beta and se is correctly computed!
MS[, `:=`(
    chr = paste("chr", CHR, sep=""),
    snp_id = paste("chr", CHR, ":", BP, sep=""),
    beta = log(OR),
    se = log(OR) / qnorm(P, lower.tail = F),
    pval = P,
    snp_ori = SNP # original
)]

MS[, `:=`(
    chr = factor(chr, autosomes),
    BEDstart = as.integer(BP) - 1,
    BEDend = as.integer(BP)
)]

# select columns, make sure to sort rows
MS_cols = c("chr", "BEDstart", "BEDend", "snp_id", "beta", "se", "pval", "snp_ori")
MS = MS[chr %in% autosomes & ! is.na(pval), ..MS_cols][order(chr, BEDstart, BEDend)]


# write 3 column bed file for coordinates
write_tsv(MS[, .(chr, BEDstart, BEDend, snp_id)],
          paste0(GWAS.out.prefix, "MS_stats_coordinates_hg19.bed"),
          col_names = F)


# Liftover from hg19 to hg38 ----------------------------------------------

#deleted.snps = fread(cmd = "awk '$1 !~ /#Dele/' /project2/yangili1/cdai/aicher/code/chRNA/Data/GWAS/MS.liftover.unmapped")

# do it in shell in the GWAS dir
# ./hg19Tohg38.sh my_hg19.bed my_hg38.bed my.unmapped


# Update hg19 cordinates in data.table to hg38 ----------------------------
MS.hg38.bed = fread(paste0(GWAS.out.prefix, "MS_stats_coordinates_hg38.bed"),
                    col.names = c("chr", "BEDstart", "BEDend", "snp_id"),
                    colClasses = c("character", rep("integer", 2), "character")) %>% unique
MS.hg38.bed = MS.hg38.bed[chr %in% autosomes] # in case some hg19 coords may be lifted to scaffolds
MS.hg38.bed[, chr := factor(chr, autosomes)]

# do inner join to bring in summary stats (data.table does right outer join)
MS38 = MS.hg38.bed[MS, on = "snp_id", nomatch = NULL][!is.na(chr)]
# remove hg19 coordinates, and old snp_id based on hg19 coords
MS38[, `:=`(i.chr = NULL, i.BEDstart = NULL, i.BEDend = NULL, snp_id = NULL)]
# reconstruct snp_id using hg38 chr and end position.
MS38[, snp_id := paste(chr, BEDend, sep=":")]
# column order, and ensure row sorted
MS38 = MS38[, ..MS_cols][order(chr, BEDstart, BEDend),] %>% unique

# now saving it to non-bed format. pos is using 1-indexed
MS38 = MS38[, .(chr, pos = BEDend, snp_id, beta, se, pval, snp_ori)]

# Finally! Write munged summary stats -------------------------------------

# write munged stats, column names in separate file
write_lines(names(MS38), paste0(GWAS.out.prefix, "MS_sumstats_hg38.colnames.txt"))
write_tsv(MS38,
          file = paste0(GWAS.out.prefix, "MS_sumstats_hg38.tsv"),
          col_names = F)



# Don't forget to bgzip and tabix in shell --------------------------------
# Also safe to delete the intermediary files: *_hg19.bed, *_hg38.bed, but keep the *.liftover.unmapped file.



