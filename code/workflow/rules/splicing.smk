'''Rules here deal with counting gene expression and splicing

Mainly: 

    - featureCounts to count gene expression
    - leafcutter to compute psi
    - deeptools to compute bamcoverage
'''

#----------------------------bigwig-------------------------------

rule BamtoBigwig:
    message: 'Make bigwig from merged Bams'
    input: 'results/{dataset}/{population}/mergedBams/{indID}_merged.bam'
    output: 'results/{dataset}/{population}/bigwig/{indID}.bw'
    log: 'logs/BamtoBigwig/{dataset}_{population}_{indID}.log'
    params: 
        TMP = "/scratch/midway2/chaodai/TMP",
        outFormat = "bigwig",
        commonFlags = " -bs 10 --effectiveGenomeSize 2862010578 --normalizeUsing RPKM ",
        filterFlags = " --minMappingQuality 20 --samFlagInclude 3 "
    threads: 8
    resources: cpu = 8, mem_mb = 25000, time = 1000
    shell: 
        '''
        bamCoverage \
            -b {input} -o {output} -of {params.outFormat} \
            -p {threads} \
            {params.commonFlags} {params.filterFlags} &> {log}
        '''


#----------------------------count expression-------------------------------

def getFeatureCountsParams(wildcards):
    if wildcards.dataset == "mRNA":
        return '-p --countReadPairs'
    elif wildcards.dataset == "chRNA":
        return '-p --countReadPairs -s 2' # same to BF's
    else: 
        sys.stderr.write("Incorrect dataset! mRNA or chRNA only.\n")
        exit(1)
    

rule featureCounts:
    input: 'results/{dataset}/{population}/mergedBams/{indID}_merged.bam'
    output: 
        summary = 'results/{dataset}/{population}/featureCounts/{indID}.counts.summary',
        counts = 'results/{dataset}/{population}/featureCounts/{indID}.counts'
    params:
        TMP = '/home/chaodai/scratch/TMP',
        prefix ='results/{dataset}/{population}/featureCounts/{indID}',
        anno = config['GENCODE_GTF'],
        options = getFeatureCountsParams
    threads: 4
    resources: cpu=4, mem_mb=15000, time=100
    shell:
        '''
        featureCounts \
            -T {threads} --tmpDir {params.TMP} \
            {params.options} \
            -a {params.anno} -o {output.counts} {input}
        ls {output.summary}
        '''

def getCollectFeatureCountsInput(wildcards):
    if wildcards.dataset == 'chRNA' and wildcards.population == 'YRI':
        infiles = expand('results/{dataset}/{population}/featureCounts/{indID}.counts',
                        dataset='chRNA', population='YRI', indID=chRNA_INDIDS)
    if wildcards.dataset == 'mRNA' and wildcards.population == 'YRI':
        infiles = expand('results/{dataset}/{population}/featureCounts/{indID}.counts',
                        dataset='mRNA', population='YRI', indID=mRNA_INDIDS)
    return infiles

rule collectFeatureCounts:
    message: 'collect feature count results into a single count table'
    input: getCollectFeatureCountsInput
    output: 'results/{dataset}/{population}/featureCounts/gathered_counts.tsv'
    log: 'logs/collectFeatureCounts/{dataset}_{population}.log'
    params:
        rscript = 'workflow/scripts/collectFeatureCounts.R'
    threads: 1
    resources: cpu = 1, mem_mb = 20000, time = 210
    shell: 
        '''
        Rscript {params.rscript} -C "{input}" -O {output} &> {log}
        '''


#----------------------------splicing-------------------------------

# I can't use merged bam because merged bam had gone through GATK's splitN
# which split all split reads into sperate reads, thus no junction exists!
# rule ExtractJunctions:
#     message: '### Extract junctions'
#     # input: 'results/{dataset}/{population}/mergedBams/{indID}_merged.bam'
#     output: 
#         junc = 'results/{dataset}/{population}/leafcutter/juncfiles/{indID}.junc'
#     params: 
#         strand = 0 # unstranded
#     resources: cpu = 1, mem = 25000, time = 1000
#     shell: 
#         '''
#         regtools junctions extract \
#             -s {params.strand} -m 50 \
#             -o {output.junc} {input}
#         '''

# rule MakeJuncFileList:
#     message: '### Make leafcutter junction files'
#     input: 
#         juncs = expand('results/leafcutter/juncfiles/{library}.junc', library=LIBRARY)
#     output: 
#         'results/leafcutter/juncfiles/juncfilelist.txt'
#     log: 'logs/MakeJuncFileList.log'
#     resources: cpu = 1, mem = 15000, time = 1000, account = 'pi-jstaley'
#     run:
#         with open(output[0], 'w') as out:
#             for filepath in input.juncs:
#                 out.write(filepath + '\n')

# rule ClusterIntrons:
#     message: '### Cluster introns'
#     input: 
#         # juncs = expand('results/leafcutter/juncfiles/{library}.junc', library=LIBRARY),
#         juncfile_list = 'results/leafcutter/juncfiles/juncfilelist.txt'
#     output: 
#         counts = 'results/leafcutter/clustering/leafcutter_perind.counts.gz',
#         numers = 'results/leafcutter/clustering/leafcutter_perind_numers.counts.gz'
#     log: 'logs/ClusterIntrons.log'
#     params: 
#         rundir = 'results/leafcutter/clustering/',
#         maxintronlength = 500000, # max intron length 500kb,
#         minclusterreads = 30 # min cluster reads
#     resources: cpu = 1, mem = 15000, time = 1000, account = 'pi-jstaley'
#     shell: 
#         '''
#         python workflow/submodules/leafcutter/scripts/leafcutter_cluster_regtools_py3.py \
#             -j {input.juncfile_list} \
#             -r {params.rundir} \
#             -l {params.maxintronlength} \
#             -m {params.minclusterreads} &> {log}
#         ls {output.counts} &> {log}
#         '''