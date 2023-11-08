'''
FURTHER FILTER A2I EDITING SITES WITHOUT REMAPPING

AUTHOR      : CHAO DAI
LAST UPDATE : 2023-09-15

This procedure takes the results from `callA2I.smk`, which has gone through 
a few rounds of basic filtering and some annotation steps, and apply further
steps of filtering. The difference in the filtering steps in this procedure
comparing to the ones in `callA2I.smk` is that these filtering steps are
generally more customized, and often require custom python or R scripts to
run. 

Like `callA2I.smk`, the lowest level for all the rules here, and onwards
are at individual/person level (`indID`).


SECTION 1
Main filtering components include:
-   Simple repeats   : Remove A2I sites that overlap with simple repeats.
-   Splice junctions : Remove A2I sites that overlap within 4 bp of splice
                       within introns.
-   Homopolyers      : Remove A2I sites that overlap within homopolymers.
                       Homopolymers are sequences with 5 or more of the same
                       bases, pre-determined from ref seq (hg38).

SECTION 2
After filtering steps, main steps include:
-   Union A2I sites : A unioned set of A2I sites from all individuals.
-   Count coverage  : Count number of reads at each A2I site.
-   Quantify EL     : Quantify (E)diting (L)evel. Here Editing level is computed
                      within each individual. The formula is:
                      EL = [ A2I reads ] / [ Total reads ] at each site. 

WILDCARDS:
-   `dataset`   : either `mRNA` or `chRNA` to indicate Geuvadis mRNA
                  dataset or in-house LCL chRNA dataset.
-   `population`: applies to both dataset, can be `YRI` or `EUR`, etc.
-   `indID`     : a person or individual ID. This is essentially the foreign key
                  used to query genotype info in 1KGP's vcf files.

NOTES ON RESULTS:
- 2023-01-14: Results of all rules in this smk file were produced before
              May 2022 using `aicher/code/mRNA-YRI` and `aichr/code/chRNA`
              snakemake rules. Results were moved into respective directory.
              Although timestamps of files from `Union_AtoI_Sites` and after
              may be refreshed `snakemake --touch`, thus reflects a more recent
              timestamp.


'''


#----------------------------------------------------------------------------
# SECTION 1: 2ND ROUND OF FILTERING OF A2I SITES
#----------------------------------------------------------------------------


# rule Filter_Primer_Error:
#     """
#     This step also infer strand for (only) A>G (T>C) mismatches.
#     """
#     message: '### Remove sites supported primer errors (6bp with 5 prim)'
#     input: 
#         TAB = "results/{dataset}/{population}/BCFCall/FinalAnno/{indID}.bed",
#         BAM = "results/{dataset}/{population}/mergedBams/{indID}_merged.bam"
#     output: 
#         tab = 'results/{dataset}/{population}/NoRemap/removePrimerError/{indID}.txt',
#         bed = 'results/{dataset}/{population}/NoRemap/removePrimerError/{indID}.bed'
#     log: 'logs/NoRemap/Filter_Primer_Error/{dataset}_{population}_{indID}.log'
#     threads: 1
#     resources: time=2000, mem_mb=35000
#     shell:
#         '''
#         python workflow/scripts/check5Preads.py \
#             --BED {input.TAB} --bam {input.BAM} --th 0.2 --outBED {output.bed} --outTAB {output.tab}  &> {log}
#         '''

rule Filter_Simple_Repeats:
    message: '### Remove simple repeats'
    input: 'results/{dataset}/{population}/BCFCall/FinalAnno/{indID}.bed'
    output: 'results/{dataset}/{population}/NoRemap/removeSimpleRepeats/{indID}.bed' 
    log: 'logs/NoRemap/Filter_Simple_Repeats/{dataset}_{population}_{indID}.log'
    params: 
        SimpleRepeats = config['SIMPLE_REPEATS']
    threads: 1
    shell:
        '''
        intersectBed -wa -v \
            -a <(cat {input} | cut -f 1-6) \
            -b {params.SimpleRepeats} 1> {output} 2>{log}
        '''


def getGatherSpliceJunctionsInput(wildcards):

    dataset, population, indID = wildcards.dataset, wildcards.population, wildcards.indID
    if dataset == 'mRNA':
        RUNIDS = mRNA_samples.query('sample_name == @indID').err_name
        SJ_FILES = expand('results/{DS}/{POP}/alignment/{RUN}.SJ.out.tab',
                           DS=dataset, POP=population, RUN=RUNIDS)
    elif dataset == 'chRNA':
        INDIDS = chRNA_samples.query('IndID == @indID').IndID
        REPS = chRNA_samples.query('IndID == @indID').RepNumber
        SJ_FILES = expand('/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/Alignments/STAR_Align/chRNA.Expression.Splicing/{IND}/{REP}/SJ.out.tab', zip, IND=INDIDS, REP=REPS)
    else:
        sys.stderr.write(f'Wildcard key error. Wildcard dataset: {dataset} is incorrect.\n')
        exit(1)
    
    return SJ_FILES

# merge SJ tab files and produce bed for later use.
rule GatherSpliceJunctions:
    message: '''### Make Bed files of coordinates of 4bp adjacent to splice sites in introns '''
    input: getGatherSpliceJunctionsInput
    output: "results/{dataset}/{population}/SJ/{indID}.SJ.4bpIntron.bed"
    log: 'logs/GatherSpliceJunctions/{dataset}_{population}_{indID}.log'
    threads: 1
    shell:
        '''
        Rscript workflow/scripts/SJ_v2.R {input} {output} &> {log}
        '''


rule Filter_Intron_SJ:
    """
    Filter intronic sites within 4bp of splice sites, 
    using specific intron splice sites for each sample produced by STAR
    """
    message: '### Remove intronic sites with 4bp of splice sites'
    input: 
        Bed = 'results/{dataset}/{population}/NoRemap/removeSimpleRepeats/{indID}.bed',
        IntronSJ = "results/{dataset}/{population}/SJ/{indID}.SJ.4bpIntron.bed"
    output: 'results/{dataset}/{population}/NoRemap/removeIntronSJ/{indID}.bed'
    log: 'logs/NoRemap/Filter_Intron_SJ/{dataset}_{population}_{indID}.log'
    threads: 1
    shell:
        '''
        intersectBed -wa -v -a {input.Bed} -b {input.IntronSJ} 1> {output} 2>{log}
        '''    


### filter homopolymers >= 5
rule Filter_homopolymers:
    message: '### Remove homopolymers of length >=5 '
    input: 'results/{dataset}/{population}/NoRemap/removeIntronSJ/{indID}.bed'
    output: 'results/{dataset}/{population}/NoRemap/removeHomopolymers/{indID}.bed'
    log: 'logs/NoRemap/Filter_Homopolymers/{dataset}_{population}_{indID}.bed'
    params: 
        Homopolymer = config['HOMOPOLYMER']
    threads: 1
    shell:
        '''
        intersectBed -wa -v -a {input} -b {params.Homopolymer} > {output}
        '''




# #----------------------------------------------------------------------------
# # SECTION 2: QUANTIFY EDITING LEVEL + EXTRA
# #----------------------------------------------------------------------------

# # Note 1: From this step onwards, editing sites are unioned across samples, thus no 
# #         longer by individual. 
# # Note 2: Previous steps include non A>I editing, but after the union step.
# #         Only A>I editing sites are kept.


def getUnion_AtoI_SitesInput(wildcards):
    dataset, population = wildcards.dataset, wildcards.population
    if dataset == 'mRNA':
        INDIDS = mRNA_INDIDS
    elif dataset == 'chRNA':
        INDIDS = chRNA_INDIDS
    else:
        sys.stderr.write(f'Wildcard key error. Wildcard dataset: {dataset} is incorrect.\n')
        exit(1)
    return expand('results/{DS}/{POP}/NoRemap/removeHomopolymers/{IND}.bed', DS=dataset, POP=population, IND=INDIDS)

rule Union_AtoI_Sites:
    """
    Output of this rule includes ONLY A>I sites (A>G, or T>C).
    """
    message: '### Union A to I sites from each individual sample'
    input: getUnion_AtoI_SitesInput
    output: 'results/{dataset}/{population}/NoRemap/UnionAISites/Unioned_A2I_sites.bed'
    log: 'logs/NoRemap/{dataset}_{population}_Union_AtoI_sites.log'
    resources: cpu = 1, mem = 15000, time = 1000
    script: '../scripts/Union_AtoI_sites.R'


rule CountCoverage:
    message: '### Count coverage at editing site'
    input: 
        bed = 'results/{dataset}/{population}/NoRemap/UnionAISites/Unioned_A2I_sites.bed',
        bam = 'results/{dataset}/{population}/mergedBams/{indID}_merged.bam'
    output: 'results/{dataset}/{population}/NoRemap/countCoverage/{indID}.txt'
    params:
        py_script = 'workflow/scripts/countCoverage_v2.py'
    log: 'logs/NoRemap/CountCoverage/{dataset}_{population}_{indID}.log'
    threads: 1
    shell:
        '''
            python {params.py_script} --BAM {input.bam} --BED {input.bed} --outTab {output} &>{log}
        '''


rule QuantEditingLevel:
    message: '### Quantify editing level'
    input: cnt = 'results/{dataset}/{population}/NoRemap/countCoverage/{indID}.txt'
    output: 'results/{dataset}/{population}/NoRemap/EditLevel/{indID}.txt'
    log: 'logs/NoRemap/QuantEditingLevel/{dataset}_{population}_{indID}.log'
    script: '../scripts/quantEditLevel.R'

rule InferStrand:
    '''Infer strand for A>G (T>C) mismatches'''
    input: rules.QuantEditingLevel.output
    output: 'results/{dataset}/{population}/NoRemap/EditLevel/{indID}_stranded.txt'
    shell:
        '''
        # use awk to add strand info to the 4th column
        # keep the first row as is
        # From second row and beyond, if the 4th column is A_G then replace 4th field as +, 
        # else replace the 4th field as -.
        # output to a new file with tab delimiter
        awk 'BEGIN {{OFS="\t"}};
             NR==1 {{print $0}};
             NR>1 {{if ($4 == "A_G") {{$4="+"}} 
                    else {{$4="-"}};
                    print $0
                  }}' {input} > {output}
        '''



def getGatherEditingInput(wildcards):
    dataset, population = wildcards.dataset, wildcards.population
    if dataset == 'mRNA':
        INDIDS = mRNA_INDIDS
    elif dataset == 'chRNA':
        INDIDS = chRNA_INDIDS
    else:
        sys.stderr.write(f'Wildcard key error. Wildcard dataset: {dataset} is incorrect.\n')
        exit(1)
    return expand('results/{DS}/{POP}/NoRemap/EditLevel/{IND}_stranded.txt', DS=dataset, POP=population, IND=INDIDS)

rule GatherEditing:
    """
    Essentially cbind all the measurement colums. 
    Basically, each PERSON file has the same rows (sites), but quantified by PERSON
    Done separately for DP: Total read counts; AP: Editing read counts; EL: editing level ([0-1], incl. NA)
    NOTE: this step DOES !!NOT!! apply any filtering.
    """
    message: '### Gather editing level, total read counts, editing read counts'
    input: getGatherEditingInput
    output: 
        DP = 'results/{dataset}/{population}/NoRemap/GatherEditing/DP.txt',
        AP = 'results/{dataset}/{population}/NoRemap/GatherEditing/AP.txt',
        EL = 'results/{dataset}/{population}/NoRemap/GatherEditing/EL.txt'
    threads: 6
    resources: cpu = 6, mem = 25000, time = 1000
    script: '../scripts/GatherEditingLevels.R'








