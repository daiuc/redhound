localrules: all

# choose these 10 libraries to subsample
LIBS = ['ERR204971', 'ERR204980', 'ERR204981', 'ERR204985', 'ERR204987', 'ERR204988', 
        'ERR204843', 'ERR204930']

CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
AUTOSOMES = [f'chr{str(i)}' for i in range(1, 23)] 

rule all:
    input:
        expand('resources/mRNA/YRI/fastq-dl/{lib}_1.fastq.gz', lib=LIBS),
        expand('resources/mRNA/YRI/fastq-dl/{lib}_2.fastq.gz', lib=LIBS)

rule SubsampleFastq:
    message:
        'Subsample fastq files to 2M reads'
    input:
        r1 = '/project2/yangili1/cdai/A2I/code/resources/mRNA/YRI/fastq-dl/{lib_name}_1.fastq.gz', 
        r2 = '/project2/yangili1/cdai/A2I/code/resources/mRNA/YRI/fastq-dl/{lib_name}_2.fastq.gz'
    output:
        r1 = 'resources/mRNA/YRI/fastq-dl/{lib_name}_1.fastq.gz',
        r2 = 'resources/mRNA/YRI/fastq-dl/{lib_name}_2.fastq.gz'
    params:
        n_reads = 500000
    threads: 1
    shell:
        '''
        seqtk sample -s100 {input.r1} {params.n_reads} | gzip -c > {output.r1}
        seqtk sample -s100 {input.r2} {params.n_reads} | gzip -c > {output.r2}
        '''

rule Extract1KGP1PercentVariants:
    message: '''Extract from 1000 Genomes Project where AF_AFR >= 0.01'''
    input: '/project2/yangili1/cdai/SNP/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz' # phase 3 1KGP
    output: temp('resources/1KGP/CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.filtered.shapeit2-duohmm-phased.1percent.vcf.gz')
    params:
        filterParams = '-i "AF_AFR >= 0.01"',
        keepSampleParams = '-s HG00096', # just select 1 sample to reduce file size
    threads: 4
    resources: cpu=4, mem_mb=16000, time=500
    shell:
        '''
        bcftools view --threads {threads} -Oz -o {output} \
            {params.filterParams} {params.keepSampleParams} \
            {input}
        '''
    

rule Concat1Extract1KGP1PercentVariants:
    message:'''Concatenate 1KGP 1% variants'''
    input: expand('resources/1KGP/CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.filtered.shapeit2-duohmm-phased.1percent.vcf.gz', chrom=AUTOSOMES)
    output: 'resources/1KGP/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1-to-22.filtered.shapeit2-duohmm-phased.1percent.vcf.gz'
    threads: 6
    resources: cpu=6, mem_mb=26000, time=1000
    shell:
        '''
        bcftools concat --threads {threads} -Oz -o {output} {input}
        bcftools index -f -t {output}
        '''


