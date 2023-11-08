library(tidyverse)
library(data.table)

runmode = 1

if (runmode == 1) {
    input_counts = snakemake@input[['cnt']]
    output = snakemake@output[[1]]
    print("snakemake script mode")
} else {
    input_counts = "Results/hs38/countCoverage/NA19130.txt"
    # input_beds = "Results/hs38/BCFCall/raw/NA18486.bed"
    output = "test.txt"
}


counts = fread(input_counts, sep = "\t", header = TRUE,
               colClasses = c(chr = "character", BEDstart = "integer",
                              BEDend = "integer", name = "character",
                              score = "integer", strand = "character",
                              A = "integer", C = "integer",
                              G = "integer", T = "integer"))

# # read bed that includes mismatch type
# beds = data.table::fread(input_beds, sep = "\t", header = FALSE)
# # keep only A>G (T>C) mismatches
# beds = filter(beds, V4 %in% c("T,C", "A,G"))

# complement counts with mismatch type
# counts = left_join(counts, beds[, 1:4],
#             by = c("chr" = "V1", "BEDstart" = "V2", "BEDend" = "V3"))  %>%
#             drop_na

# copute DP, AP, and EL (Editing Level)

# counts = rename(counts, EditType = V4)  %>%
#             mutate(DP = A + C + G + T,
#                    AP = case_when(EditType == "A,G" ~ G,
#                                   EditType == "T,C" ~ C))  %>% 
#             mutate(EL = round(AP/DP, 5))  %>%
#             select(-A, -C, -G, -T) %>%
#             mutate_at(c("BEDstart", "BEDend"), as.integer)

counts[, `:=`(DP = A + C + G + T,
              AP = if_else(name == "A_G", G, C)
)]
counts[, EL := round(AP/DP, 5)]

# write to tsv
counts[, .(chr, BEDstart, BEDend, name, score, strand, DP, AP, EL)]  %>%
    fwrite(., output, sep = "\t", na = "NA")