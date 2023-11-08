library(tidyverse)
library(data.table)


# set run mode: 1 - snakemake ; 2 - debug in rstudio directly
runmode = 1

# inputfile is tab delimited output in
if (runmode == 1) {
    input.files = snakemake@input
    output.file = snakemake@output[[1]]
    print("Running in snakemake script mode")
} else {
    input.files = dir('Results/hs38/removeHomopolymers', "^[A-Z]{2}[0-9]+\\.bed$", full.names = T)
}

# read in bed files
names(input.files) = str_extract(input.files, "[A-Z]{2}[0-9]{4,6}")
print("Union sites from these samples:")
print(names(input.files))
beds = map(input.files, ~ fread(.x))
beds = rbindlist(beds)

# chroms
CHROMS = c(paste("chr", 1:22, sep=""), "chrX", "chrY")

# select on A>G or T>C sites
beds[, V1 := factor(V1, CHROMS)]
beds = beds[V4 %in% c("A_G", "T_C")][, .(V1, V2, V3, V4, V5=999, V6)] %>%
        unique %>% 
        .[order(V1, V2)]

fwrite(beds, output.file, sep = "\t", col.names = FALSE)
