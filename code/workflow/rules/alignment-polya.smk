'''
PolyA-RNA-seq starts from fastq files

Everything runs at library or run ID level. Each individual (person) can have
multiple `lib`. 

wildcards:
    - population: population name, etc. YRI
    - lib: library or run ID, etc. ERR188021

'''


#----------------------------------------------------------------------------
# POLYA-RNA-SEQ SEPCIFIC RULES
#----------------------------------------------------------------------------

# Note atomic level of wildcards here is a sequencing run. A sample can have 
# multiple runs, and each run can have read1 and read2 fastqs.

rule RemoveAdapters:
    message: '''### Trim adapters ###'''
    input:
        f1 = "resources/mRNA/{population}/fastq-dl/{lib}_1.fastq.gz",
        f2 = "resources/mRNA/{population}/fastq-dl/{lib}_2.fastq.gz"
    output:
        f1 = temp("resources/mRNA/{population}/trimmed-fastq/{lib}_R1.fastq.gz"),
        f2 = temp("resources/mRNA/{population}/trimmed-fastq/{lib}_R2.fastq.gz"),
        report_html = "resources/mRNA/{population}/trimmed-fastq/{lib}.fastp_report.html",
        report_json = "resources/mRNA/{population}/trimmed-fastq/{lib}.fastp_report.json"
    log: 'logs/RemoveAdapter/mRNA_{population}_{lib}.log'
    wildcard_constraints:
        lib = 'ERR\d+'
    threads: 1
    resources: time=200, mem_mb=20000, cpu=1
    shell:
        '''
        fastp --in1 {input.f1} --out1 {output.f1} \
            --in2 {input.f2} --out2 {output.f2} \
            --thread {threads} \
            --dont_overwrite \
            --overrepresentation_analysis \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --disable_quality_filtering \
            --length_required 15 \
            --json {output.report_json} \
            --html {output.report_html} &> {log}
        '''


# # STAR alignment, output sorted bam
rule STAR:
    message: '''### Alignment using STARS ###'''
    input:
        f1 = rules.RemoveAdapters.output.f1,
        f2 = rules.RemoveAdapters.output.f2
    output: 
        bam = "results/mRNA/{population}/alignment/{lib}.Aligned.sortedByCoord.out.bam",
        SJ = 'results/mRNA/{population}/alignment/{lib}.SJ.out.tab'
    wildcard_constraints:
        lib = 'ERR\d+'
    params:
        genome_dir = config['STAR_INDEX_DIR'],
        out_prefix = 'results/mRNA/{population}/alignment/{lib}.',
        read_files_comm = "zcat",
        VCF = config['STAR_VCF'],
        limitOutSJcollapsed = 5000000, # adding this for nuclear RNA, otherwise report error sometimes
        outSAMattrRGline = "ID:{lib} SM:{lib} PL:ILLUMINA LB:{lib} PU:{lib}",
        outSAMattributes = "NH HI AS NM MD RG vW",
        outBAMparams = "--outFileNamePrefix results/mRNA/{population}/alignment/{lib}. --outSAMtype BAM SortedByCoordinate --outSAMunmapped None"
    threads: 6
    resources: time=2100, mem_mb=40000, cpu=6
    log: 'logs/STAR/mRNA_{population}_{lib}.log'
    #     BAM SortedByCoordinate output Aligned.sortedByCoord.out.bam file
    #     BAM Unsorted output Aligned.out.bam
    shell:
        '''
        STAR --runThreadN {threads} --genomeDir {params.genome_dir} --twopassMode Basic \
            --readFilesIn {input.f1} {input.f2} --readFilesCommand {params.read_files_comm} \
            {params.outBAMparams} \
            --waspOutputMode SAMtag --varVCFfile <(zcat {params.VCF}) \
            --outSAMattrRGline {params.outSAMattrRGline} --outSAMattributes {params.outSAMattributes} \
            --outFilterType BySJout --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.1 \
            --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
            --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
            --limitOutSJcollapsed {params.limitOutSJcollapsed} &> {log}
        
        samtools index -@ {threads} {output.bam} &>> {log}
 
        '''

rule FilterBam:
    message: '''### Filter bam file ###'''
    input:
        bam = rules.STAR.output.bam
    output:
        bam = temp("results/mRNA/{population}/alignment/{lib}.Filtered.bam")
    params:
        regions = ' '.join(CHROMS),
        expr = '-e " ![vW] || [vW] == 1 "'
    log: 'logs/FilterBam/mRNA_{population}_{lib}.log'
    threads: 4
    resources: time=200, mem_mb=25000, cpu=4
    shell:
        '''
        samtools view -@ {threads} -b -F 4 -F 256 {params.expr} \
            -o {output.bam} {input.bam} {params.regions} &> {log}
        '''


rule MarkDups:
    message: '''### MarkDuplications - mRNA ###'''
    input: 
        bam = rules.FilterBam.output.bam
    output:
        bam = temp("results/mRNA/{population}/alignment/{lib}.Sorted.Deduped.bam"),
        metrics = "results/mRNA/{population}/alignment/{lib}.Dedup.Metrics"
    params:
        tmp = "/scratch/midway2/chaodai/TMP" # added because gatk's tmp directory error
    log: 'logs/MarkDups/mRNA_{population}_{lib}.log'
    threads: 1
    resources: time=2000, mem_mb=20000, cpu=1
    shell:
        '''
        gatk MarkDuplicates  \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --TMP_DIR {params.tmp} \
            --REMOVE_SEQUENCING_DUPLICATES &> {log}
        '''


# # BQSR on RNA-seq requires first split reads that are mapped to multiple loci, 
# # indicated by the N number of cigar string. 
rule SplitNCigarReads:
    message: '''### Split chimeric RNA-seq reads using the N Cigar string ###'''
    input:
        bam = rules.MarkDups.output.bam,
        ref = config['FA_HS38']
    output:
        bam = temp("results/mRNA/{population}/BQSR/{lib}.splitNCigar.bam")
    params:
        tmp = "/scratch/midway2/chaodai/TMP" # added because gatk's tmp directory error
    log: 'logs/SplitNCigarReads/mRNA_{population}_{lib}.log'
    threads: 1
    resources: time=2000, mem_mb=20000, cpu=1
    shell:
        '''
        gatk SplitNCigarReads \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.bam} \
            --tmp-dir {params.tmp} &> {log}
            
            # --skip-mapping-quality-transform by default, gatk will change mapping quality to 60
        '''


# BQSR
# Note GATK BaseRecalibrator will not take in any reads that has N in cigar
rule BaseRecalibrator:
    message: '''### Compute covariate matrix of Base Recalibration ###'''
    input:
        bam = rules.SplitNCigarReads.output.bam,
    output:
        bam = "results/mRNA/{population}/BQSR/{lib}_recal.bam",
        recal_file = "results/mRNA/{population}/BQSR/{lib}_covariates.tab",
        flag = touch("results/mRNA/{population}/BQSR/{lib}_recal.done")
    params:
        ref = config['FA_HS38'],
        known_sites = config['dbSNP'],
        tmp = "/scratch/midway2/chaodai/TMP" # added because gatk's tmp directory errorthreads: 1
    log: 'logs/BaseRecalibrator/mRNA_{population}_{lib}.log'
    resources: time=2000, mem_mb=20000, cpu=1
    shell:
        '''
        gatk BaseRecalibrator -I {input.bam} -R {params.ref} --known-sites {params.known_sites} \
            --tmp-dir {params.tmp} \
            -O {output.recal_file} &> {log}
        
        sleep 5 
        
        gatk ApplyBQSR \
            -I {input.bam} \
            -R {params.ref} \
            --bqsr-recal-file {output.recal_file} \
            --tmp-dir {params.tmp} \
            -O {output.bam} &>> {log}

        '''



