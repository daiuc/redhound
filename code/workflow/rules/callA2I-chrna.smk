'''
Call RNA editing sites for chRNA-seq data

'''

# call variants
use rule CallVariant as CallVariant_chrna with:
    message: '''### Call variants on chRNA###'''
    input:
        bam = "results/chRNA/{population}/BQSR/{indID}.addRG.Deduped.splitNCigar.recal.bam"
    output: 
        vcf = "results/chRNA/{population}/BCFCall/raw/{indID}.vcf.gz"
    params: 
        ref = config['FA_HS38'],
        region = "" # use "-r chr21" for limited regions
    log: 'logs/CallVariant/chRNA_{population}_{indID}.log'
    threads: 4
    resources: time=2000, mem_mb=20000, cpu=4


# filter variants using these criteria: 
# 1. SNPs only (remove indels)
# DP (total read depth) > 9, 
# within chr1 - chrY
# QUAL > 19 
use rule RawFilterVCF as RawFilterVCF_chrna with:
    input: 
        vcf = rules.CallVariant_chrna.output.vcf
    output: "results/chRNA/{population}/BCFCall/rawfilters/{indID}.vcf.gz"
    params: 
        regions = ",".join(CHROMS),
        include_expr = "'INDEL=0 && QUAL > 19 && DP>9 && INFO/AD[1]>0'"
    log: 'logs/RawFilterVCF/chRNA_{population}_{indID}.log'


use rule RemoveCommonVariants as RemoveCommonVariants_chrna with:
    input: rules.RawFilterVCF_chrna.output
    output: "results/chRNA/{population}/BCFCall/remove_common/{indID}.vcf.gz"
    params:
        tmp = '/scratch/midway3/chaodai/TMP/chRNA_{population}_BCFCall_RemoveVariants_{indID}',
        KGP3_1PERCENT = config['KGP3_AF_AFR_GT_1PCT'] # 1KGP3 variants AF_AFR > 1%
    log: 'logs/RemoveCommonVariants/chRNA_{population}_{indID}.log'

use rule TabulateVariants as TabulateVariants_chrna with:
    input: rules.RemoveCommonVariants_chrna.output
    output: "results/chRNA/{population}/BCFCall/remove_common/{indID}.tsv.gz"
    params:
        fields = lambda f: "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/DP\t%INFO/AD{0}\t%INFO/AD{1}\n", 
        header = '"CHROM\tPOS\tID\tREF\tALT\tDP\tAD_ref\tAD_alt"'
    log: 'logs/TabulateVariants/chRNA_{population}_{indID}.log'




