if (exists("snakemake")) {
    INPUTS = snakemake@input
    OUT.DP = snakemake@output[['DP']]
    OUT.AP = snakemake@output[['AP']]
    OUT.EL = snakemake@output[['EL']]
    NCores = snakemake@threads[[1]]
    print("Working with snakemake script")
} else {
    INPUTS = scan(text = "Results/hs38/EditLevel/NA19150.txt Results/hs38/EditLevel/NA19200.txt",
                  what = "character")
    NCores = 8
}

NCores = as.numeric(NCores)
library(furrr)
plan(multisession, workers = min(8, availableCores(), NCores))
options(future.globals.maxSize = 26843545600) # set to 25G
print(paste0("Using multisession with furrr with ", min(8, availableCores(), NCores), " threads"))

library(tidyverse)
library(data.table)


names(INPUTS) = str_extract(INPUTS, "[A-Z]{2}[0-9]{5}")

dt = future_map(INPUTS,
    ~ fread(.x, sep = "\t", header = TRUE,
            colClasses = c("character", "integer", "integer",
                           "character", "integer", "character",
                           "integer", "integer", "double")))

# check rows are exactly in the same order across samples
row_nms = future_map(dt, 
    ~ .x[, .(rn = paste(chr, BEDstart, BEDend, name, sep = "_"))] %>%
        .$rn)

rowCheck = future_map_lgl(row_nms, ~ all(row_nms[[1]] == .x))
rowCheck = all(rowCheck)

if (rowCheck) {
    print("All rows match. Proceed...")
    if (exists("row_nms")) rm(row_nms) # collect memory

    SharedDT = dt[[1]][, .(chr, BEDstart, BEDend, name, score, strand)]
    dt.DP = future_map(dt, ~ .x[, DP]) %>% do.call(cbind, .) %>% as.data.table
    dt.AP = future_map(dt, ~ .x[, AP]) %>% do.call(cbind, .) %>% as.data.table
    dt.EL = future_map(dt, ~ .x[, EL]) %>% do.call(cbind, .) %>% as.data.table

    dt.DP = cbind(SharedDT, dt.DP)
    dt.AP = cbind(SharedDT, dt.AP)
    dt.EL = cbind(SharedDT, dt.EL)

    walk2(list(dt.DP, dt.AP, dt.EL),
          list(OUT.DP, OUT.AP, OUT.EL),
          ~ write_tsv(.x, .y))

    print(paste0("Finished at: ", date()))
} else {
    print("Error. Something happened. Go figure.")
}
