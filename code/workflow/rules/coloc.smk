'''
COLOCALIZATIONS

AUTHOR      : CHAO DAI
LAST UPDATE : 2023-01-14

This procedure performs colocalization between edQTLs and GWAS traits. Further colocalization
analyses of other molecular traits, such as eQTLs, spQTLs with GWAS, or with edQTLs should
also be appended is smk file.

GWAS SUMMARY STATISTICS:
-   The original GWAS summary dataset are saved at `resources/GWAS/gwas_download`.
    Cleaned up GWAS summary dataset are saved at `resources/GWAS/{gwas}_sumstats_hg38.tsv.gz`.
-   Column numbers in each sumstats can differ, but same columns (e.g. 'chr', 'pos') will
    have consistent column names.
-   Summary stats are compressed with `bgzip` and indexed with `tabix` for fast retrieval.
-   Descriptions of the GWAS studies should be updated in `resources/GWAS/gwas_config.yaml`.


SECTION 1: COLOC - EDQTLS VS. GWAS

WILDCARDS:
-   `dataset`   : either `mRNA` or `chRNA` to indicate Geuvadis mRNA
                  dataset or in-house LCL chRNA dataset.
-   `population`: applies to both dataset, can be `YRI` or `EUR`, etc.
-   `gwas`      : GWAS trait name, e.g. 'Height'.



'''




#----------------------------------------------------------------------------
# SECTION 1: COLOC - EDQTLS VS. GWAS
#----------------------------------------------------------------------------


rule IdentifyGwasHitLoci:
    '''
    Current bugs: 
        Overlapping hit loci are merged into one, using `reduce`, which can 
        result in large (multi Mbp) locus. This is too large!
    
    '''
    message: '### Find GWAS hit loci from sumstats. These are used as coloc focus windows.'
    input: 
        sumstats = 'resources/GWAS/{gwas}_sumstats_hg38.tsv.gz',
        colnames = 'resources/GWAS/{gwas}_sumstats_hg38.colnames.txt'
    output: 'resources/GWAS/{gwas}_Hits_1e-7_hg38.tsv'
    log: 'logs/IdentifyGwasHitLoci/{gwas}.log'
    params: 
        pval_th = 1e-7 # pvalue threshold
    threads: 1
    resources: cpu = 1, mem_mb = 20000, time = 2100
    script: '../scripts/GwasHits-v2.R'


rule GWAS_v_edQTL:
    message: '### Run colocolization analysis'
    input: 
        GWAS_sumstats = 'resources/GWAS/{gwas}_sumstats_hg38.tsv.gz',
        GWAS_sumstats_cols = 'resources/GWAS/{gwas}_sumstats_hg38.colnames.txt',
        GWAS_hits = 'resources/GWAS/{gwas}_Hits_1e-7_hg38.tsv', # GWAS hits use 1Mbp windows
        GWAS_config = 'resources/GWAS/gwas_config.yaml',
        QTL_sumstats = 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_1000000/{chrom}/Nominal.txt.gz', # QTL sumstats use 1Mbp cis window
        QTL_hits = 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_10000/permutations_pass_all_chr.AddedQvalue.txt.gz' # QTL hits use 10kb cis window
    output: touch('results/{dataset}/{population}/coloc/{gwas}/{chrom}/coloc.done')
    params: 
        coloc_summary = 'results/{dataset}/{population}/coloc/{gwas}/{chrom}/summary.tsv',
        QTL_threshold = 0.1,
        earlystop = 'results/{dataset}/{population}/coloc/{gwas}/{chrom}/StoppedEarly.txt',
    threads: 1
    resources: cpu = 1, mem = 30000, time = 2000
    script: '../scripts/coloc_V2.R'

rule GatherColocSummary:
    message: 'Gather coloc summary'
    input: 
        GWAS_sumstats = 'resources/GWAS/{gwas}_sumstats_hg38.tsv.gz',
        GWAS_sumstats_cols = 'resources/GWAS/{gwas}_sumstats_hg38.colnames.txt',
        dbsnp150_lookup = "/project2/yangili1/cdai/SNP/snp150_hg38_3cols.txt.gz"
    output: 'results/{dataset}/{population}/coloc/{gwas}/coloc_loci_summary.tsv'
    params: 
        QTL_sumstats_prefix = 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_1000000/',
        COLOC_summary_prefix = 'results/{dataset}/{population}/coloc/{gwas}/chr'
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2100
    script: "../scripts/collectColocResults.R"


rule ColocWithEqtl:
    '''
    Make sure to check all inputs and params, particularly column names must match
    '''
    message: '### Run coloc: edQTL & eQTL'
    wildcard_constraints: # restrictive for now
        dataset = 'chRNA',
        population = 'YRI'
    input: 
        perm1 = 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_10000/permutations_pass_all_chr.AddedQvalue.txt.gz',
        perm2 = 'resources/bjf79-qtltools/chRNA.Expression.Splicing/PermutationPassForColoc.txt.gz'
    output: 'results/{dataset}/{population}/coloc/eQTL/coloc_summary.tsv'
    params:
        perm1cols = 'pid,pchr,pstart,pend,pstrand,n_snps,best_dist,best_gid,best_gchr,best_gstart,best_gend,dof_true,dof_est,beta_ml1,beta_ml2,pval_nom,pval_r2,slope,pval_emp,pval_beta_adj,q',
        perm2cols = 'pid,pchr,pstart,pend,pstrand,n_snps,best_dist,best_gid,best_gchr,best_gstart,best_gend,dof_true,dof_est,beta_ml1,beta_ml2,pval_nom,pval_r2,slope,slope_se,pval_emp,pval_beta_adj',
        nom1fileprefix = 'results/chRNA/YRI/edQTL_FixedPCs/ranknorm/cis_1000000/', # prefix may be file
        nom2fileprefix = 'resources/bjf79-qtltools/chRNA.Expression.Splicing/NominalPassForColoc.txt.tabix.gz', # prefix may be file
        nom1cols = 'pid,pchr,pstart,pend,pstrand,n_snps,dist,snp_id,gchr,gstart,gend,pval.nom,pval.r2,beta,topSNP',
        nom2cols = 'pid,pchr,pstart,pend,pstrand,n_snps,dist,snp_id,gchr,gstart,gend,pval.nom,pval.r2,beta,beta_se,topSNP',
        rscript = 'workflow/scripts/coloc_molQTLs.R'
    threads: 8
    resources: cpu = 8, mem_mb = 25000, time = 2100
    log: 'logs/ColocWithEqtl/{dataset}_{population}.log'
    shell: 
        '''
        Rscript {params.rscript} -c {threads} \
            --perm1 {input.perm1} --perm1cols {params.perm1cols} \
            --nom1 {params.nom1fileprefix} --nom1cols {params.nom1cols} \
            --perm2 {input.perm2} --perm2cols {params.perm2cols} \
            --nom2 {params.nom2fileprefix} --nom2cols {params.nom2cols} \
            -O {output} &> {log}
        '''

use rule ColocWithEqtl as ColocWithSqtl with:
    message: '### Run coloc: edQTL & sQTL'
    input:
        perm1 = 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_10000/permutations_pass_all_chr.AddedQvalue.txt.gz',
        perm2 = 'resources/bjf79-qtltools/chRNA.Splicing/PermutationPassForColoc.txt.gz'
    output: 'results/{dataset}/{population}/coloc/sQTL/coloc_summary.tsv'
    params:
        perm1cols = 'pid,pchr,pstart,pend,pstrand,n_snps,best_dist,best_gid,best_gchr,best_gstart,best_gend,dof_true,dof_est,beta_ml1,beta_ml2,pval_nom,pval_r2,slope,pval_emp,pval_beta_adj,q',
        perm2cols = 'pid,pchr,pstart,pend,pstrand,n_snps,best_dist,best_gid,best_gchr,best_gstart,best_gend,dof_true,dof_est,beta_ml1,beta_ml2,pval_nom,pval_r2,slope,slope_se,pval_emp,pval_beta_adj',
        nom1fileprefix = 'results/chRNA/YRI/edQTL_FixedPCs/ranknorm/cis_1000000/', # prefix may be file
        nom2fileprefix = 'resources/bjf79-qtltools/chRNA.Splicing/NominalPassForColoc.txt.tabix.gz', # prefix may be file
        nom1cols = 'pid,pchr,pstart,pend,pstrand,n_snps,dist,snp_id,gchr,gstart,gend,pval.nom,pval.r2,beta,topSNP',
        nom2cols = 'pid,pchr,pstart,pend,pstrand,n_snps,dist,snp_id,gchr,gstart,gend,pval.nom,pval.r2,beta,beta_se,topSNP',
        rscript = 'workflow/scripts/coloc_molQTLs.R'
    log: 'logs/ColocWithSqtl/{dataset}_{population}.log'

use rule ColocWithEqtl as ColocWithNcqtl with:
    message: '### Run coloc: edQTL & ncRNA QTL'
    input:
        perm1 = 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_10000/permutations_pass_all_chr.AddedQvalue.txt.gz',
        perm2 = 'resources/bjf79-qtltools/chRNA.Expression_ncRNA/PermutationPassForColoc.txt.gz'
    output: 'results/{dataset}/{population}/coloc/ncQTL/coloc_summary.tsv'
    params:
        perm1cols = 'pid,pchr,pstart,pend,pstrand,n_snps,best_dist,best_gid,best_gchr,best_gstart,best_gend,dof_true,dof_est,beta_ml1,beta_ml2,pval_nom,pval_r2,slope,pval_emp,pval_beta_adj,q',
        perm2cols = 'pid,pchr,pstart,pend,pstrand,n_snps,best_dist,best_gid,best_gchr,best_gstart,best_gend,dof_true,dof_est,beta_ml1,beta_ml2,pval_nom,pval_r2,slope,slope_se,pval_emp,pval_beta_adj',
        nom1fileprefix = 'results/chRNA/YRI/edQTL_FixedPCs/ranknorm/cis_1000000/', # prefix may be file
        nom2fileprefix = 'resources/bjf79-qtltools/chRNA.Expression_ncRNA/NominalPassForColoc.txt.tabix.gz', # prefix may be file
        nom1cols = 'pid,pchr,pstart,pend,pstrand,n_snps,dist,snp_id,gchr,gstart,gend,pval.nom,pval.r2,beta,topSNP',
        nom2cols = 'pid,pchr,pstart,pend,pstrand,n_snps,dist,snp_id,gchr,gstart,gend,pval.nom,pval.r2,beta,beta_se,topSNP',
        rscript = 'workflow/scripts/coloc_molQTLs.R'
    log: 'logs/ColocWithNcqtl/{dataset}_{population}.log'

rule AnnotateColoc:
    message: 'Annotate colocalization summary results'
    input: 
        coloc = 'results/{dataset}/{population}/coloc/{qtltype}/coloc_summary.tsv',
        siteAnno = 'results/{dataset}/{population}/GatherEditing/AnnotatedSites.bed',
        ncRNA = config['NCRNA']
    output: 'results/{dataset}/{population}/coloc/{qtltype}/coloc_summary_annotated.tsv'
    wildcard_constraints: # restrictive for now
        dataset = 'chRNA',
        population = 'YRI',
        qtltype = 'eQTL|sQTL|ncQTL'
    log: 'logs/Annotatecoloc/{dataset}_{population}_{qtltype}.log'
    params:
        annoCols = "chr,start,end,edit,score,strand,genefeature,genename,genetype,repeat,db",
        rscript = 'workflow/scripts/annoColocOutput.R'
    threads: 8
    resources: cpu = 8, mem_mb = 20000, time = 2100
    shell: 
        '''
        Rscript {params.rscript} -c {threads} \
            -I {input.coloc} \
            -E {input.siteAnno} --annoCols {params.annoCols} \
            -N {input.ncRNA} \
            -O {output} &> {log}
        '''       


def getLocusCompareDataParams(wildcards):
    nom2fileReplace = {
        'eQTL': 'chRNA.Expression.Splicing',
        'sQTL': 'chRNA.Splicing',
        'ncQTL': 'chRNA.Expression_ncRNA'
    }
    nom2fileprefix = f'resources/bjf79-qtltools/{nom2fileReplace[wildcards.qtltype]}/NominalPassForColoc.txt.tabix.gz'
    return nom2fileprefix

rule getLocusCompareData:
    input: 'results/{dataset}/{population}/coloc/{qtltype}/coloc_summary_annotated.tsv',
    output: 
        json = 'results/{dataset}/{population}/coloc/{qtltype}/locusCompareData.json'
    wildcard_constraints: # restrictive for now
        dataset = 'chRNA',
        population = 'YRI',
        qtltype = 'eQTL|sQTL|ncQTL'
    threads: 1
    resources: cpu = 1, mem_mb = 30000, time = 1000
    log: 'logs/getLocusCompareData/{dataset}_{population}_{qtltype}.log'
    params:
        nom1fileprefix = 'results/chRNA/YRI/edQTL_FixedPCs/ranknorm/cis_1000000/', # prefix may be file
        nom2fileprefix = getLocusCompareDataParams, # prefix may be file
        nom1cols = 'pid,pchr,pstart,pend,pstrand,n_snps,dist,snp_id,gchr,gstart,gend,pval.nom,pval.r2,beta,topSNP',
        nom2cols = 'pid,pchr,pstart,pend,pstrand,n_snps,dist,snp_id,gchr,gstart,gend,pval.nom,pval.r2,beta,beta_se,topSNP',
        minH4PP = 0.7,
        rscript = 'workflow/scripts/getLocusCompareData.R'
    threads: 1
    resources: cpu = 1, mem_mb = 20000, time = 2100
    shell:
        '''
           Rscript {params.rscript} -c {threads} \
                --colocFile {input} \
                --nom1 {params.nom1fileprefix} \
                --nom1cols {params.nom1cols} \
                --nom2 {params.nom2fileprefix} \
                --nom2cols {params.nom2cols} \
                --minH4PP {params.minH4PP} \
                -O {output} &> {log}
        '''
    

def getMakeLocusComparePlotsPerm2(wildcards):
    folder = {
    'eQTL': 'chRNA.Expression.Splicing',
    'sQTL': 'chRNA.Splicing',
    'ncQTL': 'chRNA.Expression_ncRNA'
    }
    perm2 = f'resources/bjf79-qtltools/{folder[wildcards.qtltype]}/PermutationPassForColoc.txt.gz'
    return perm2

rule makeLocusComparePlots:
    input: 
        json = 'results/{dataset}/{population}/coloc/{qtltype}/locusCompareData.json',
        coloc = 'results/{dataset}/{population}/coloc/{qtltype}/coloc_summary_annotated.tsv',
        perm1 = 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_10000/permutations_pass_all_chr.AddedQvalue.txt.gz',
        perm2 = getMakeLocusComparePlotsPerm2
    output: 'results/{dataset}/{population}/coloc/{qtltype}/locusComparePlots.pdf'
    wildcard_constraints: # restrictive for now
        dataset = 'chRNA',
        population = 'YRI',
        qtltype = 'eQTL|sQTL|ncQTL'
    params:
        perm1cols = 'pid,pchr,pstart,pend,pstrand,n_snps,best_dist,best_gid,best_gchr,best_gstart,best_gend,dof_true,dof_est,beta_ml1,beta_ml2,pval_nom,pval_r2,slope,pval_emp,pval_beta_adj,q',
        perm2cols = 'pid,pchr,pstart,pend,pstrand,n_snps,best_dist,best_gid,best_gchr,best_gstart,best_gend,dof_true,dof_est,beta_ml1,beta_ml2,pval_nom,pval_r2,slope,slope_se,pval_emp,pval_beta_adj',
        rscript = 'workflow/scripts/locusComparePlots.R'
    threads: 1
    resources: cpu = 1, mem_mb = 30000, time = 1000
    log: 'logs/makeLocusComparePlots/{dataset}_{population}_{qtltype}.log'
    shell:
        '''
        Rscript {params.rscript} \
            --jsonFile {input.json} --colocFile {input.coloc} \
            --perm1 {input.perm1} --perm1cols "{params.perm1cols}" \
            --perm2 {input.perm2} --perm2cols "{params.perm2cols}" \
            -O {output} &> {log}
        '''




## pygenome tracks
rule getPlotRegionFile:
    input:
        coloc_sum = 'results/{dataset}/{population}/coloc/{qtltype}/coloc_summary_annotated.tsv',
        plot_regions = 'resources/230220-edQTL-{qtltype}-LocusComp-Handpicked-ids.txt'
    output: 'results/{dataset}/{population}/coloc/{qtltype}/pygenometrack_regions.tsv'
    params:
        rscript = 'workflow/scripts/getGenomeTrackRegions.R'
    shell:
        '''
        Rscript {params.rscript} -I {input.coloc_sum} -R {input.plot_regions} -O {output}
        '''



