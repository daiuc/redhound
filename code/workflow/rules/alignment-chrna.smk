'''
LCL CHRNA-seq starts from bam files from Ben's directory

Every rule here runs at the individual (person) level.

wildcards:
    - population: population name, etc. YRI
    - indID: individual ID, etc. NA19130

'''


def getMergeBamByIndInput(wildcards):
    indID = wildcards.indID
    REPS = chRNA_samples.query('IndID == @indID').RepNumber.unique()
    return expand('test-chrna-bams/{indID}/{Rep}/Filtered.bam', indID=indID, Rep=REPS)
    # return expand("/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/Alignments/STAR_Align/chRNA.Expression.Splicing/{indID}/{Rep}/Filtered.bam", indID=indID, Rep=REPS)


rule mergeBamByInd:
    message: '### merge bam files by individual'
    input: getMergeBamByIndInput
    output: 
        bam = temp('results/chRNA/{population}/alignment/{indID}.bam'),
        bai = temp('results/chRNA/{population}/alignment/{indID}.bam.bai')
    log: 'logs/mergeBamByInd/chRNA_{population}_{indID}.log'
    threads: 6
    resources: cpu=6, mem_mb=25000, time=500
    shell:
        '''
        samtools merge -@ {threads} - {input} 2> {log} | samtools sort -@ {threads} -o {output.bam} - 2>> {log}
        samtools index {output.bam} &>> {log}
        ls {output} &>> {log}
        '''


rule addReadGroups:
    message: '### Add RG tags to chRNA bams'
    input:
        bam = rules.mergeBamByInd.output.bam
    output: 
        bam = temp("results/chRNA/{population}/alignment/{indID}.addRG.bam"),
        bai = temp("results/chRNA/{population}/alignment/{indID}.addRG.bam.bai")
    log: "logs/addReadGroups/chRNA_{population}_{indID}.log"
    params:
        RGID = '{indID}',
        RGSM = '{indID}',
        RGPL = "ILLUMINA",
        RGLB = '{indID}',
        RGPU = '{indID}' 
    threads: 2
    resources: time = 500, mem_mb = 30000, cpu = 2
    shell:
        '''
        gatk AddOrReplaceReadGroups -I {input.bam} -O {output.bam} \
            --RGID {params.RGID} --RGSM {params.RGSM} --RGPL {params.RGPL} --RGLB {params.RGLB} --RGPU {params.RGPU} &> {log}
        samtools index -@ {threads} {output.bam} &>> {log}
        ls {output} &>> {log}
        '''


use rule MarkDups as MarkDups_chrna with:
    message: '### Mark duplicates in chRNA-seq bam files'
    input:
        bam = rules.addReadGroups.output.bam
    output:
        bam = temp("results/chRNA/{population}/alignment/{indID}.addRG.Deduped.bam"),
        metrics = "results/chRNA/{population}/alignment/{indID}.addRG.Deduped.Metrics"
    params:
        tmp = "/scratch/midway2/chaodai/TMP" # added because gatk's tmp directory error
    log: "logs/MarkDups/chRNA_{population}_{indID}.log"
    threads: 1
    resources: time=2000, mem_mb=30000, cpu=1


use rule SplitNCigarReads as SplitNCigarReads_chrna with:
    message: '### Split chimeric RNA-seq reads using the N Cigar string'
    input:
        bam = rules.MarkDups_chrna.output.bam,
        ref = config['FA_HS38']
    output:
        bam = temp("results/chRNA/{population}/BQSR/{indID}.addRG.Deduped.splitNCigar.bam"),
        bai = temp("results/chRNA/{population}/BQSR/{indID}.addRG.Deduped.splitNCigar.bai")
    params:
        tmp = "/scratch/midway2/chaodai/TMP" # added because gatk's tmp directory error
    log: "logs/SplitNCigarReads/chRNA_{population}_{indID}.log"
    threads: 1
    resources: time=2000, mem_mb=30000, cpu=1



use rule BaseRecalibrator as BaseRecalibrator_chrna with:
    message: '### BQSR for bam files ###'
    input:
        bam = rules.SplitNCigarReads_chrna.output.bam,
    output:
        bam = "results/chRNA/{population}/BQSR/{indID}.addRG.Deduped.splitNCigar.recal.bam",
        recal_file = "results/chRNA/{population}/BQSR/{indID}_covariates.tab",
        flag = touch("results/chRNA/{population}/BQSR/{indID}_recal.done")
    params:
        ref = config['FA_HS38'],
        known_sites = config['dbSNP'],
        tmp = "/scratch/midway2/chaodai/TMP" # added because gatk's tmp directory errorthreads: 1
    log: "logs/BaseRecalibrator/chRNA_{population}_{indID}.log"
    resources: time=2000, mem_mb=30000, cpu=1





# #----------------------------------------------------------------------------
# # SECTION 2: LCL CHRNA SEPCIFIC RULES
# #----------------------------------------------------------------------------
