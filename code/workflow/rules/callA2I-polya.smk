'''
Call RNA editing sites from polyA RNA-seq data


wildcards:
    - population: population name, etc. YRI
    - indID: individual ID, etc. NA19130, can corresponds to multiple libraries
    - lib: library or run ID, etc. ERR188021

'''

def getMergeRecaledBamInputs(wildcards):
    LIBS = mRNA_samples.query('sample_name == @wildcards.indID').err_name.unique()
    LIBS = [lib for lib in LIBS if lib in test_m_libs] # comment out when not testing
    return expand("results/mRNA/{pops}/BQSR/{lib}_recal.bam", pops=wildcards.population, lib=LIBS)


rule MergeRecaledBam:
    message: '''### Merge and filter splitted and BQSR applied bam files per individual ###'''
    input: 
        bams = getMergeRecaledBamInputs,
    output: 
        bam = "results/mRNA/{population}/BQSR/{indID}.Filtered.Merged.Recaled.bam"
    params:
        samtoolsMergeParams = "-f --write-index", 
        samtoolsViewParams = "-h -b --write-index",
        tmp = '/scratch/midway3/chaodai/TMP',
        RGID = '{indID}',
        RGSM = '{indID}',
        RGPL = "ILLUMINA",
        RGLB = '{indID}',
        RGPU = '{indID}' 
    log: 'logs/MergeRecaledBam/mRNA_{population}_{indID}.log'
    threads: 4
    resources: time=500, mem_mb=20000, cpu=4
    shell:
        '''
        tmp_bam1={params.tmp}/{wildcards.indID}.tmp1.bam

        samtools merge -@ {threads} {params.samtoolsMergeParams} \
            -o $tmp_bam1 {input.bams} &> {log} # merge

        gatk AddOrReplaceReadGroups -I $tmp_bam1 -O {output.bam} \
            --RGID {params.RGID} --RGSM {params.RGSM} --RGPL {params.RGPL} --RGLB {params.RGLB} --RGPU {params.RGPU} &> {log}
        
        samtools index -@ {threads} {output.bam} &>> {log} 

        rm {params.tmp}/{wildcards.indID}.tmp*
        '''



# call variants
rule CallVariant:
    message: '''### Call variants ###'''
    input:
        bam = rules.MergeRecaledBam.output.bam
    output: 
        vcf = "results/mRNA/{population}/BCFCall/raw/{indID}.vcf.gz"
    params: 
        ref = config['FA_HS38'],
        region = "" # use "-r chr21" for limited regions
    log: 'logs/CallVariant/mRNA_{population}_{indID}.log'
    threads: 4
    resources: time=2000, mem_mb=20000, cpu=4
    shell:
        '''
        bcftools mpileup --threads {threads} \
            -q 20 -Q 20 {params.region} --skip-any-set DUP \
            -a FORMAT/AD,FORMAT/SP,INFO/AD \
            -Ou -f {params.ref} {input.bam} 2> {log} | \
            bcftools call -m -v --group-samples - -Oz -o {output.vcf} --threads {threads} &>> {log}
            bcftools index -f -t {output.vcf} &>> {log}
        '''

# filter variants using these criteria: 
# 1. SNPs only (remove indels)
# DP (total read depth) > 9, 
# within chr1 - chrY, this is taken care of at mergeRecalBam step
# QUAL > 19 
rule RawFilterVCF:
    message: '''### Apply 1st set of filters ###'''
    input:
        vcf = rules.CallVariant.output.vcf
    output: "results/mRNA/{population}/BCFCall/rawfilters/{indID}.vcf.gz"
    params:
        regions = ",".join(CHROMS),
        include_expr = "'INDEL=0 && QUAL > 19 && DP>9 && INFO/AD[1]>0'"
    log: 'logs/RawFilterVCF/mRNA_{population}_{indID}.log'
    threads: 2
    resources: time=2000, mem_mb=10000, cpu=2
    shell:
        '''
        bcftools filter \
            -i {params.include_expr} -r {params.regions}\
            --threads {threads}  \
            -Oz -o {output} {input.vcf} &> {log}
        bcftools index -f -t {output} &>> {log}
        
        '''

rule RemoveCommonVariants:
    message: '''### Remove variants that overlap with common genetic variants ###'''
    input: rules.RawFilterVCF.output
    output: "results/mRNA/{population}/BCFCall/remove_common/{indID}.vcf.gz"
    params:
        tmp = '/scratch/midway3/chaodai/TMP/mRNA_{population}_BCFCall_RemoveVariants_{indID}',
        KGP3_1PERCENT = config['KGP3_AF_AFR_GT_1PCT'] # 1KGP3 variants AF_AFR > 1%
    log: 'logs/RemoveCommonVariants/mRNA_{population}_{indID}.log'
    threads: 2
    resources: time=2000, mem_mb=15000, cpu=2
    shell:
        '''
        mkdir -p {params.tmp}
        bcftools isec -p {params.tmp} -Oz -n~10 {input} {params.KGP3_1PERCENT} &> {log}
        mv {params.tmp}/0000.vcf.gz {output}
        bcftools index -f -t {output} &>> {log}
        rm -rf {params.tmp}
        '''


rule TabulateVariants:
    message: '''### Output variants into tabulated text file ###'''
    input: rules.RemoveCommonVariants.output
    output: "results/mRNA/{population}/BCFCall/remove_common/{indID}.tsv.gz"
    log: 'logs/TabulateVariants/mRNA_{population}_{indID}.log'
    params:
        fields = lambda f: "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/DP\t%INFO/AD{0}\t%INFO/AD{1}\n", 
        header = '"CHROM\tPOS\tID\tREF\tALT\tDP\tAD_ref\tAD_alt"'
    shell:
        '''
        bcftools query -f '{params.fields}' {input} 2> {log} | \
            awk '
                BEGIN {{OFS="\t"; print {params.header}}};
                      {{print $0}}
                ' 2>> {log} | bgzip -c 1> {output} 2>> {log}
        '''




#---------------------------------------------------------------
#   Alternatively, use HaplotypeCaller to call variants
#---------------------------------------------------------------

rule HaplotypeCaller:
    input: bam = rules.MergeRecaledBam.output.bam
    output: 
        gvcf = "results/mRNA/{population}/HaplotypeCaller/raw/{indID}.g.vcf.gz"
    params: ref = config['FA_HS38']
    log: 'logs/HaplotypeCaller/mRNA_{population}_{indID}.log'
    shell:
        '''
        gatk HaplotypeCaller \
            -R {params.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ERC GVCF &> {log}
        '''

rule GenotypeGVCFs:
    input: rules.HaplotypeCaller.output.gvcf
    output: vcf = "results/mRNA/{population}/HaplotypeCaller/raw/{indID}.vcf.gz"
    log: 'logs/GenotypeGVCFs/mRNA_{population}_{indID}.log'
    params: 
        ref = config['FA_HS38'],
        tmp = '/scratch/midway3/chaodai/TMP'
    shell:
        '''
        gatk GenotypeGVCFs \
            -R {params.ref} \
            -V {input} \
            -O {output.vcf} \
            --tmp-dir {params.tmp} &> {log}
        '''

