library(tidyverse)
library(data.table)
library(RNOmni)

if (exists("snakemake")) {
    print("Using snakemake:")
    print(str(snakemake))

    DP = snakemake@input[['DP']]
    AP = snakemake@input[['AP']]
    EL = snakemake@input[['EL']]
    ANNO = snakemake@input[['ANNO']]

    SAMPLE_OUT = snakemake@output[[1]]
    AP_raw_OUT = snakemake@output[[2]]
    AP_norm_OUT = snakemake@output[[3]]
    EL_raw_OUT = snakemake@output[[4]]
    EL_norm_OUT = snakemake@output[[5]]

    MIN_N_SAMPLES = as.numeric(snakemake@params[['MIN_N_SAMPLES']])
    MIN_DP_PER_SAMPLE = as.numeric(snakemake@params[['MIN_DP_PER_SAMPLE']])
    MIN_AP_PER_SAMPLE = as.numeric(snakemake@params[['MIN_AP_PER_SAMPLE']])
    MIN_AP_ROWSUM = as.numeric(snakemake@params[['MIN_AP_ROWSUM']])
    EL_ROWMEAN_RANGE = scan(text = snakemake@params[['EL_ROWMEAN_RANGE']], what="character") %>% as.numeric
    EXCLUDE_SAMPLES = snakemake@params[['EXCLUDE_SAMPLES']]

} else {
    INPUT_FILES = scan(text = "Results/hs38/GatherEditing/DP.txt Results/hs38/GatherEditing/AP.txt Results/hs38/GatherEditing/EL.txt",
                 what = "character")
    DP = INPUT_FILES[[1]]
    AP = INPUT_FILES[[2]]
    EL = INPUT_FILES[[3]]
    MIN_N_SAMPLES = 10
    MIN_DP_PER_SAMPLE = 10
    MIN_AP_PER_SAMPLE = 1
    MIN_AP_ROWSUM = 25
    EL_ROWMEAN_RANGE = c(1e-3, 0.99)
    EXCLUDE_SAMPLES = "NA18855"
    ANNO = 'Results/hs38/GatherEditing/AnnotatedSites.bed'
}


# Compute the number of samples passing given threshold for a given site
compute_N_per_th = function(m, th) {
    #' @param m: DP or AP matrix
    #' @param th: thresholds,a numeric vector
    m = rowSums(modify(m, ~ .x >= th))
    return(m)
}

# imput NAs
fill_NAs = function(m) {
    #' @param m: EL matrix
    m = apply(m, 1, function(row) {
        row[which(is.na(row))] = mean(row, na.rm = T) %>% round(5)
        return(row)
    })
    return(t(m))
}

myscale <- function(v) {
    # handle vectors with SD=0 
  if (sd(v, na.rm = T) == 0) {
    v = c(v, mean(v, na.rm = T) + 1e-5)
    v_scaled = scale(v)[, 1]
    v_scaled = v_scaled[1:length(v_scaled)-1]
  } else {
    v_scaled = scale(v)[, 1]
  }
  
  return(v_scaled)
}

##################
print("### Load in DP, AP, and EL")

headers = read_lines(DP, n_max = 1)  %>% str_split("\t")  %>% unlist
coltypes_DP_AP = c("character", "integer", "integer", "character", "integer", "character",
             rep("integer", length(headers) - 6))
coltypes_EL = c("character", "integer", "integer", "character", "integer", "character",
             rep("double", length(headers) - 6))

DP = fread(DP, sep = "\t", header = TRUE, colClasses = coltypes_DP_AP)
AP = fread(AP, sep = "\t", header = TRUE, colClasses = coltypes_DP_AP)
EL = fread(EL, sep = "\t", header = TRUE, colClasses = coltypes_EL)


# select sites (rows) each having at least MIN_N_SAMPLES satisfying DP and AP requirements
print(paste0("### Subset sites using Min_N_Samples satisfying requirements, DP: ",
        MIN_DP_PER_SAMPLE, " , AP: ", MIN_AP_PER_SAMPLE,
        " , >= ", MIN_N_SAMPLES, " samples"))

DP.idx = compute_N_per_th(DP[, -c(1:6)], th = MIN_DP_PER_SAMPLE)  %>% `>=`(MIN_N_SAMPLES) %>% which
AP.idx = compute_N_per_th(AP[, -c(1:6)], th = MIN_AP_PER_SAMPLE)  %>% `>=`(MIN_N_SAMPLES) %>% which
idx = intersect(DP.idx, AP.idx)
DP = DP[idx,]
AP = AP[idx,]
EL = EL[idx,]
rm(idx)


# select sites(rows) each passing rowsum(AP) and  rowmean(EL) requirement
print(paste0("### Subset sites using rowsums of AP ",  MIN_AP_PER_SAMPLE,
     " and rowmeans of EL: ", EL_ROWMEAN_RANGE[1], " , ", EL_ROWMEAN_RANGE[2]))

rowsum.idx = which(rowSums(AP[, -c(1:6)]) >= MIN_AP_ROWSUM)
rowmean.idx = which(rowMeans(EL[, -c(1:6)], na.rm = T) > EL_ROWMEAN_RANGE[1] &
               rowMeans(EL[, -c(1:6)], na.rm = T) < EL_ROWMEAN_RANGE[2])
idx = intersect(rowsum.idx, rowmean.idx)
DP = DP[idx,]
AP = AP[idx,]
EL = EL[idx,]


# prepare QTLtools phenotype data [BED]
SAMPLE_NAMES = colnames(DP)[str_detect(colnames(DP), "[A-Z]{2}[0-9]{5}")]
print(paste0("### Removing sample ", EXCLUDE_SAMPLES))
SAMPLE_NAMES = SAMPLE_NAMES[! SAMPLE_NAMES %in% EXCLUDE_SAMPLES]
print("Remaining samples:")
print(SAMPLE_NAMES)

CHROMS = c(paste("chr", 1:22, sep = ""), "chrX", "chrY")

shared_dt = DP[, .(chr, BEDstart, BEDend, name, score, strand)]
shared_dt[, `:=`(
    pid = paste(chr, BEDend, strand, name, sep = "_")
)]

# load annotations
print("### Load annotated sites to get gid")
# print(ANNO)
ANNO = fread(ANNO, sep = "\t", header = F,
             col.names = c("chr", "BEDstart", "BEDend", "name", "score", "strand",
                           "feature", "gene", "gene_type", "repeat", "radar"))
# ensures protein coding genes ranked first
ANNO[, `:=`(
    gene_type = factor(gene_type, c("protein_coding", "lncRNA", "miRNA", "snRNA", "snoRNA", ".")),
    feature = factor(feature, c("UTR", "exon", "intron", "."))
)]

# perform ROJ to bring in annotation, this introduces duplicates of rows
# because some pid may have multiple genes
# just rank and pick one
lookup = ANNO[, -c("name", "score")][
            shared_dt, on = c("chr", "BEDstart", "BEDend", "strand")]

lookup[, `:=`(
    rk = rank(gene_type, feature, ties.method = "first"),
    n = .N
), by = .(pid)]

# remove duplicates
lookup = lookup[rk == 1, -c('rk', 'n')]
shared_dt = lookup[, -c("chr", "BEDstart", "BEDend", "strand", "name", "score")][
    shared_dt, on = "pid"]

if (all(shared_dt$BEDstart == EL$BEDstart)) {
    print("All common rows match. Outputing phenotype bed files.")

    # first row normalization to N(0,1)
    # AP.mx.norm = AP[, ..SAMPLE_NAMES] %>% t %>% scale %>% t
    # EL.mx.norm = fill_NAs(EL[, ..SAMPLE_NAMES]) %>% t %>% scale %>% t

    # using custom scale function to handle vectors with SD=0
    AP.mx.norm = AP[, ..SAMPLE_NAMES] %>% t %>% as.data.frame %>% map_dfc(~myscale(.x)) %>% t
    EL.mx.norm = fill_NAs(EL[, ..SAMPLE_NAMES]) %>% t %>% as.data.frame %>% map_dfc(~myscale(.x)) %>% t

    # then Inverse Normal Transform columns and round to 5 digits
    AP.mx.norm = apply(AP.mx.norm, 2, RankNorm) %>% apply(., 2, round, 5)
    EL.mx.norm = apply(EL.mx.norm, 2, RankNorm) %>% apply(., 2, round, 5)

    # convert matrix to data.table
    AP.mx.norm = as.data.table(AP.mx.norm)
    EL.mx.norm = as.data.table(EL.mx.norm)
    colnames(AP.mx.norm) = SAMPLE_NAMES
    colnames(EL.mx.norm) = SAMPLE_NAMES

    # prep data.table
    pheno_AP_bed.raw = bind_cols(shared_dt[, .(chr, BEDstart, BEDend, pid, gene, strand)],
              AP[, ..SAMPLE_NAMES])
    pheno_AP_bed.norm = bind_cols(shared_dt[, .(chr, BEDstart, BEDend, pid, gene, strand)],
                                  AP.mx.norm)
    setorder(pheno_AP_bed.raw, chr, BEDstart, BEDend, strand)
    setorder(pheno_AP_bed.norm, chr, BEDstart, BEDend, strand)
    names(pheno_AP_bed.raw) = c("#Chr", "start", "end", "pid", "gid", "strand", SAMPLE_NAMES)
    names(pheno_AP_bed.norm) = c("#Chr", "start", "end", "pid", "gid", "strand", SAMPLE_NAMES)

    pheno_EL_bed.raw = bind_cols(shared_dt[, .(chr, BEDstart, BEDend, pid, gene, strand)],
                             EL[, ..SAMPLE_NAMES])
    pheno_EL_bed.norm = bind_cols(shared_dt[, .(chr, BEDstart, BEDend, pid, gene, strand)],
                                  EL.mx.norm)
    setorder(pheno_EL_bed.raw, chr, BEDstart, BEDend, strand)
    setorder(pheno_EL_bed.norm, chr, BEDstart, BEDend, strand)
    names(pheno_EL_bed.raw) = c("#Chr", "start", "end", "pid", "gid", "strand", SAMPLE_NAMES)
    names(pheno_EL_bed.norm) = c("#Chr", "start", "end", "pid", "gid", "strand", SAMPLE_NAMES)

    # write out
    write_tsv(pheno_AP_bed.raw, AP_raw_OUT)
    write_tsv(pheno_AP_bed.norm, AP_norm_OUT)
    write_tsv(pheno_EL_bed.raw, EL_raw_OUT)
    write_tsv(pheno_EL_bed.norm, EL_norm_OUT)
    write_lines(SAMPLE_NAMES, file = SAMPLE_OUT)

    print("Finished.")
}
