library(data.table)
library(tidyverse)



if(interactive()){
    print("Running in interactive mode")
    args = "/project2/yangili1/cdai/aicher/code/chRNA/Results/hs38/QTLprep/scale_center/EL_raw.bed.gz"
} else{
    print("Running in snakemake script mode")
    args <- commandArgs(trailingOnly=TRUE)
}


FILEIN = args[[1]]
FILEOUT = args[[2]]

dt = fread(FILEIN)
old_col_names = names(dt)

shuffled_sample_cols = sample(7:ncol(dt))
while (any(shuffled_sample_cols == 7:ncol(dt))) {
    # continue shuffling until no column is in the same order
    shuffled_sample_cols = sample(7:ncol(dt))
}

dt = dt[, c(1:6, shuffled_sample_cols), with=F]
names(dt) = old_col_names

write_tsv(dt, FILEOUT)
