'''

'''


#----------------------------------------------------------------------------
# SECTION 1: 2ND ROUND OF FILTERING OF A2I SITES
#----------------------------------------------------------------------------


def getMergeVCFsInput(wildcards):
    '''Returns all A2I editing sites of 1 population of 1 dataset
    '''
    population = wildcards.population
    INDIDS = ['NA19095', 'NA19185'] # modify for real data!
    return expand('results/mRNA/{POP}/BCFCall/remove_common/{IND}.vcf.gz', POP=population, IND=INDIDS)

rule MergeVCFs:
    input: getMergeVCFsInput
    output: 'results/mRNA/{population}/BCFCall/merged.vcf.gz'
    log: 'logs/MergeVCFs/mRNA_{population}.log'
    threads: 4
    resources: time=500, mem_mb=30000, cpu=4
    shell:
        '''
        bcftools merge --threads {threads} -Oz -o {output} {input} &> {log}
        bcftools index --threads {threads} -t {output} &>> {log}
        '''




# rule Gather1stPassSites:
#     """
#     Gather all the editing sites from first pass (result of `callA2I.smk`).
#     This results in a VCF file that's compatable for STAR wasp.
#     """
#     message: '''### Make VCF file from A>G editing sites'''
#     input: getGather1stPassSitesInput
#     output: "results/mRNA/{population}/Remap/STAR-remap-variant.vcf"
#     threads: 1
#     resources: time=120, mem_mb=8000, cpu=1
#     script: "scripts/makeA2G_VCF.R"


# def Allocate_Bam_to_Fastq_Mem(wildcards, extra = 1000):
#     """
#     Allocate appropriate memory based on file size
#     Memory size in unit of MB
#     """
#     import os
#     file = "Results/hs38/mergedBams/" + wildcards.indID + "_merged.bam"
#     fileSize = os.path.getsize(file)
#     mem_mb = int(round(fileSize / 1e9) * 1000 + extra)
#     return mem_mb

# This step requires a lot of tmp space when bam files are large, e.g. with chRNA
# rule ExtractFastqFromBam:
#     message: '''### Extract paired end fastq from BAM'''
#     input: 'results/mRNA/{population}/mergedBams/{indID}_merged.bam'
#     output: 
#         R1 = temp("results/mRNA/{population}/BamToFastq/{indID}.R1.b2f.fastq.gz"),
#         R2 = temp("results/mRNA/{population}/BamToFastq/{indID}.R2.b2f.fastq.gz"),
#     log: 'logs/ExtractFastqFromBam/mRNA_{population}_{indID}.log'
#     params:
#         collate_tmp = "/scratch/midway2/chaodai/TMP/mRNA/{population}/{indID}.collate.tmp.bam",
#         prefix = "/scratch/midway2/chaodai/TMP/mRNA/{population}/{indID}"
#     threads: 4
#     resources: 
#         time=2000, 
#         mem_mb=16000, 
#         cpu=4
#     shell:
#         '''
#         samtools collate -@ {threads} -o {params.collate_tmp} {input} {params.prefix} &> {log} 
#         samtools fastq -@ {threads} -1 {output.R1} -2 {output.R2} -0 /dev/null -s /dev/null {params.collate_tmp} &>> {log}
#         rm {params.collate_tmp}
#         '''


# rule Remap_STAR:
#     """
#     Remap reads that overlap 1st pass A>G sites
#     """
#     message: '''### Remap using STARS ###'''
#     input:
#         R1 = rules.ExtractFastqFromBam.output.R1,
#         R2 = rules.ExtractFastqFromBam.output.R2,
#         A2G_VCF = "results/mRNA/{population}/Remap/STAR-remap-variant.vcf"
#     output: temp("results/mRNA/{population}/Remap/{indID}.bam")
#     log: 'logs/Remap_STAR/mRNA_{population}_{indID}.log'
#     params:
#         GENOMEDIR = config['STAR_INDEX_DIR'],
#         REF_FA = config['FA_HS38'],
#         OUT_PREFIX = "results/mRNA/{population}/Remap/{indID}.",
#         read_files_comm = "zcat",
#         limitOutSJcollapsed = 5000000, # adding this for nuclear RNA, otherwise report error sometimes
#     threads: 12
#     resources: time=2000, mem_mb=38000, cpu=12, partition="broadwl", account="pi-jstaley"
#     shell:
#         '''
#         STAR="/software/STAR-2.7.7a-el7-x86_64/bin/STAR"
#         $STAR \
#             --runThreadN {threads} \
#             --genomeDir {params.GENOMEDIR} \
#             --readFilesIn {input.R1} {input.R2} \
#             --readFilesCommand {params.read_files_comm} \
#             --outFilterType BySJout \
#             --outSAMattributes NH HI AS NM MD RG vW \
#             --outSAMunmapped None \
#             --outFilterMultimapNmax 20 \
#             --alignSJoverhangMin 8 \
#             --alignSJDBoverhangMin 1 \
#             --outFilterMismatchNmax 999 \
#             --outFilterMismatchNoverReadLmax 0.1 \
#             --alignIntronMin 20 \
#             --alignIntronMax 1000000 \
#             --alignMatesGapMax 1000000 \
#             --outSAMtype BAM Unsorted \
#             --outFileNamePrefix {params.OUT_PREFIX} \
#             --limitOutSJcollapsed {params.limitOutSJcollapsed} \
#             --outSAMattrRGline ID:{wildcards.indID} SM:{wildcards.indID} PL:ILLUMINA LB:{wildcards.indID} PU:{wildcards.indID} \
#             --varVCFfile {input.A2G_VCF} \
#             --waspOutputMode SAMtag &> {log}
        
#         sleep 10
#         samtools sort -@ {threads} -o {output} {params.OUT_PREFIX}Aligned.out.bam &>> {log}
#         sleep 10
#         samtools index -@ {threads} {output} &>> {log}
#         sleep 10
#         rm {params.OUT_PREFIX}Aligned.out.bam &>> {log}
        
#         '''


# ### filter remapped reads to include only reads that pass WASP
# rule Filter_Remapped_Bam:
#     message: '### Remove alignments failing to pass WASP tag in the 2nd STAR mapping'
#     input: "results/mRNA/{population}/Remap/{indID}.bam"
#     output: "results/mRNA/{population}/Remap/filterBam/{indID}.bam"
#     log: 'logs/Filter_Remapped_Bam/results/mRNA_{population}_{indID}.log'
#     threads: 4
#     resources: time=2000, mem_mb=20000, cpu=4, partition="broadwl"
#     shell:
#         '''
#         samtools view -@ {threads} -h -F 256 -F 2048 -q 20 \
#           -e " [vW] == 1" \
#           -o {output} {input} &> {log}
#         sleep 10
#         samtools index -@ {threads} {output} &>> {log}
#         '''

# ### filter editing sites based on remap consistency 
# rule Filter_Remap_Consistency:
#     message: '### Filter sites based on remap requirements'
#     input:
#         BED = "results/mRNA/{population}/BCFCall/FinalAnno/{indID}.bed",
#         BAMPRE = 'results/mRNA/{population}/mergedBams/{indID}_merged.bam',
#         BAMPOST = "results/mRNA/{population}/Remap/filterBam/{indID}.bam"
#     output:
#         tab = "results/mRNA/{population}/RemapConsistency/{indID}.txt",
#         bed = "results/mRNA/{population}/RemapConsistency/{indID}.bed"
#     log: 'logs/Remap_Consistency/mRNA_{population}_{indID}.log'
#     threads: 1
#     resources: time=2000, mem_mb=35000, cpu=1, partition="broadwl"
#     shell:
#         '''
#         python workflow/scripts/checkRemap.py \
#             --BED {input.BED} --bamPre {input.BAMPRE} --bamPost {input.BAMPOST} \
#             --outTab {output.tab} --outBED {output.bed} &> {log}
#         '''

# rule Filter_Primer_Error:
#     """
#     This step also infer strand for (only) A>G (T>C) mismatches.
#     """
#     message: '### Remove sites supported primer errors (6bp with 5 prim)'
#     input: 
#         TAB = "results/mRNA/{population}/RemapConsistency/{indID}.txt",
#         BAM = "results/mRNA/{population}/Remap/filterBam/{indID}.bam"
#     output: 
#         tab = 'results/mRNA/{population}/removePrimerError/{indID}.txt',
#         bed = 'results/mRNA/{population}/removePrimerError/{indID}.bed'
#     log: 'logs/Filter_Primer_Error/mRNA_{population}_{indID}.log'
#     threads: 1
#     resources: time=2000, mem_mb=35000
#     shell:
#         '''
#         python workflow/scripts/check5Preads.py \
#             --BED {input.TAB} --bam {input.BAM} --th 0.2 --outBED {output.bed} --outTAB {output.tab}  &> {log}
#         '''

# rule Filter_Simple_Repeats:
#     message: '### Remove simple repeats'
#     input: 'results/mRNA/{population}/removePrimerError/{indID}.bed'
#     output: 'results/mRNA/{population}/removeSimpleRepeats/{indID}.bed' 
#     log: 'logs/Filter_Simple_Repeats/mRNA_{population}_{indID}.log'
#     params: 
#         SimpleRepeats = config['SIMPLE_REPEATS']
#     threads: 1
#     shell:
#         '''
#         intersectBed -wa -v -a {input} -b {params.SimpleRepeats} 1> {output} 2>{log}
#         '''

# def getGatherSpliceJunctionsInput(wildcards):

#     dataset, population, indID = wildcards.dataset, wildcards.population, wildcards.indID
#     if dataset == 'mRNA':
#         RUNIDS = mRNA_samples.query('sample_name == @indID').err_name
#         SJ_FILES = expand('results/{DS}/{POP}/alignment/{RUN}.SJ.out.tab',
#                            DS=dataset, POP=population, RUN=RUNIDS)
#     elif dataset == 'chRNA':
#         INDIDS = chRNA_samples.query('IndID == @indID').IndID
#         REPS = chRNA_samples.query('IndID == @indID').RepNumber
#         SJ_FILES = expand('/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/Alignments/STAR_Align/chRNA.Expression.Splicing/{IND}/{REP}/SJ.out.tab', zip, IND=INDIDS, REP=REPS)
#     else:
#         sys.stderr.write(f'Wildcard key error. Wildcard dataset: mRNA is incorrect.\n')
#         exit(1)
    
#     return SJ_FILES

# # merge SJ tab files and produce bed for later use.
# rule GatherSpliceJunctions:
#     message: '''### Make Bed files of coordinates of 4bp adjacent to splice sites in introns '''
#     input: getGatherSpliceJunctionsInput
#     output: "results/mRNA/{population}/SJ/{indID}.SJ.4bpIntron.bed"
#     log: 'logs/GatherSpliceJunctions/mRNA_{population}_{indID}.log'
#     threads: 1
#     shell:
#         '''
#         Rscript workflow/scripts/SJ_v2.R {input} {output} &> {log}
#         '''


# rule Filter_Intron_SJ:
#     """
#     Filter intronic sites within 4bp of splice sites, 
#     using specific intron splice sites for each sample produced by STAR
#     """
#     message: '### Remove intronic sites with 4bp of splice sites'
#     input: 
#         Bed = 'results/mRNA/{population}/removeSimpleRepeats/{indID}.bed',
#         IntronSJ = "results/mRNA/{population}/SJ/{indID}.SJ.4bpIntron.bed"
#     output: 'results/mRNA/{population}/removeIntronSJ/{indID}.bed'
#     log: 'logs/Filter_Intron_SJ/mRNA_{population}_{indID}.log'
#     threads: 1
#     shell:
#         '''
#         intersectBed -wa -v -a {input.Bed} -b {input.IntronSJ} 1> {output} 2>{log}
#         '''    


# ### filter homopolymers >= 5
# rule Filter_homopolymers:
#     message: '### Remove homopolymers of length >=5 '
#     input: 'results/mRNA/{population}/removeIntronSJ/{indID}.bed'
#     output: 'results/mRNA/{population}/removeHomopolymers/{indID}.bed'
#     log: 'logs/Filter_Homopolymers/mRNA_{population}_{indID}.bed'
#     params: 
#         Homopolymer = config['HOMOPOLYMER']
#     threads: 1
#     shell:
#         '''
#         intersectBed -wa -v -a {input} -b {params.Homopolymer} > {output}
#         '''




# #----------------------------------------------------------------------------
# # SECTION 2: QUANTIFY EDITING LEVEL + EXTRA
# #----------------------------------------------------------------------------

# # Note 1: From this step onwards, editing sites are unioned across samples, thus no 
# #         longer by individual. 
# # Note 2: Previous steps include non A>I editing, but after the union step.
# #         Only A>I editing sites are kept.


# def getUnion_AtoI_SitesInput(wildcards):
#     dataset, population = wildcards.dataset, wildcards.population
#     if dataset == 'mRNA':
#         INDIDS = mRNA_INDIDS
#     elif dataset == 'chRNA':
#         INDIDS = chRNA_INDIDS
#     else:
#         sys.stderr.write(f'Wildcard key error. Wildcard dataset: mRNA is incorrect.\n')
#         exit(1)
#     return expand('results/{DS}/{POP}/removeHomopolymers/{IND}.bed', DS=dataset, POP=population, IND=INDIDS)

# rule Union_AtoI_Sites:
#     """
#     Output of this rule includes ONLY A>I sites (A>G, or T>C).
#     """
#     message: '### Union A to I sites from each individual sample'
#     input: getUnion_AtoI_SitesInput
#     output: 'results/mRNA/{population}/UnionAISites/Unioned_A2I_sites.bed'
#     log: 'logs/mRNA_{population}_Union_AtoI_sites.log'
#     resources: cpu = 1, mem = 15000, time = 1000
#     script: '../scripts/Union_AtoI_sites.R'


# rule CountCoverage:
#     message: '### Count coverage at editing site'
#     input: 
#         bed = 'results/mRNA/{population}/UnionAISites/Unioned_A2I_sites.bed',
#         bam = 'results/mRNA/{population}/Remap/filterBam/{indID}.bam'
#     output: 'results/mRNA/{population}/countCoverage/{indID}.txt'
#     log: 'logs/CountCoverage/mRNA_{population}_{indID}.log'
#     threads: 1
#     shell:
#         '''
#             python workflow/Scripts/countCoverage_v2.py \
#                 --BAM {input.bam} --BED {input.bed} --outTab {output} &>{log}
#         '''


# rule QuantEditingLevel:
#     message: '### Quantify editing level'
#     input: cnt = 'results/mRNA/{population}/countCoverage/{indID}.txt'
#     output: 'results/mRNA/{population}/EditLevel/{indID}.txt'
#     log: 'logs/QuantEditingLevel/mRNA_{population}_{indID}.log'
#     script: '../scripts/quantEditLevel.R'


# def getGatherEditingInput(wildcards):
#     dataset, population = wildcards.dataset, wildcards.population
#     if dataset == 'mRNA':
#         INDIDS = mRNA_INDIDS
#     elif dataset == 'chRNA':
#         INDIDS = chRNA_INDIDS
#     else:
#         sys.stderr.write(f'Wildcard key error. Wildcard dataset: mRNA is incorrect.\n')
#         exit(1)
#     return expand('results/{DS}/{POP}/EditLevel/{IND}.txt', DS=dataset, POP=population, IND=INDIDS)

# rule GatherEditing:
#     """
#     Essentially cbind all the measurement colums. 
#     Basically, each PERSON file has the same rows (sites), but quantified by PERSON
#     Done separately for DP: Total read counts; AP: Editing read counts; EL: editing level ([0-1], incl. NA)
#     NOTE: this step DOES !!NOT!! apply any filtering.
#     """
#     message: '### Gather editing level, total read counts, editing read counts'
#     input: getGatherEditingInput
#     output: 
#         DP = 'results/mRNA/{population}/GatherEditing/DP.txt',
#         AP = 'results/mRNA/{population}/GatherEditing/AP.txt',
#         EL = 'results/mRNA/{population}/GatherEditing/EL.txt'
#     threads: 6
#     resources: cpu = 6, mem = 25000, time = 1000
#     script: '../scripts/GatherEditingLevels.R'

# rule AnnotateEditing:
#     """
#     Annotate editing sites with gene names, introns, repfamily
#     columns are: chrom start end edit score strand feature gene_name gene_type repeat_family db
#     """
#     message: '### Annotate editing sites with gene features, gene names, RADAR db, repeats'
#     input: 'results/mRNA/{population}/GatherEditing/DP.txt'
#     output: 'results/mRNA/{population}/GatherEditing/AnnotatedSites.bed'
#     params:
#         Exons_Introns_Utrs = config['EXONS_INTRONS_UTRS'],
#         Repeats = config['REPEATS'],
#         Radar = config['DB']
#     resources: cpu = 1, mem = 15000, time = 600
#     log: 'logs/mRNA_{population}_AnnotateEditing.log'
#     shell:
#         '''
#         (awk 'BEGIN {{ OFS="\t" }}; NR > 1 {{ print $1,$2,$3,$4,$5,$6 }}' {input} | \
#             intersectBed -loj -s -a stdin -b {params.Exons_Introns_Utrs} | cut -f 1-6,10,13,14 | \
#             intersectBed -loj -s -a stdin -b {params.Repeats} | cut -f 1-9,13 | \
#             intersectBed -loj -a stdin -b {params.Radar} | cut -f 1-10,14 > {output} ) &> {log}
#         '''






