'''
CALL EDITING QTLS

AUTHOR      : CHAO DAI
LAST UPDATE : 2023-01-14

This procedure follows `filterA2I.smk`. It prepares all the necessary phenotype, genotype, and
covariance matrices for QTL calling. QTLs are called using qtltools.


SECTION 1: PREP PHENOTYPE, GENOTYPE, COVARIANCE MATRIX
-   Phenotype preparation
-   Genotype preparation
-   Prepare covariance matrix


# SECTION 2: MAP QTLS (PERMUTATION PASS)
-   Map edQTLs using qtltools permutation pass by chromosome, chr1 -> chr22
-   Combine results from each chromosome
-   Compute and add FDR (qvalue)

SECTION 3: MAP QTLS (PERMUTATION) USING RANDOMIZED PHENOTYPE AS CONTROL
-   Randomize phenotype bed file
-   Map edQTLs on randomized phenotype

# SECTION 4: MAP QTLs (nominal pass) 

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
              Although timestamps of files in this smk file may be refreshed
              `snakemake --touch`, thus reflects a more recent timestamp.
              
              Previous `aicher` runs tried using qtltools' own `--normal` flag
              to normalize editing level. These QTLs results were saved under
              `edQTL_FixedPCs/scale_center` directory. Both nominal and
              permutation pass were run in the `aicher` pipeline. Results were
              still kept in the `aicher/code/chRNA` directory.

              Permuted phenotype based QTL calls were only run against the
              chRNA dataset, not the mRNA dataset. Because of this, when writing
              the `all` rule, do not include outputs of these extra rules.




'''

#----------------------------------------------------------------------------
# SECTION 1: PREP PHENOTYPE, GENOTYPE, COVARIANCE MATRIX
#----------------------------------------------------------------------------

rule prep_phenotype:
    '''
    directory `QTLprep/ranknorm` for ranknorm normalized editing level
    directory `QTLrpep/scale_center` for non-normalized raw editing level,
    which is later normalized in qtltools using --normal flag.
    Notice a number of filters used in `params`.
    '''
    message: '### Preparing Phenotype BED file and Sample list'
    input: 
        DP = 'results/{dataset}/{population}/GatherEditing/DP.txt',
        AP = 'results/{dataset}/{population}/GatherEditing/AP.txt',
        EL = 'results/{dataset}/{population}/GatherEditing/EL.txt',
        ANNO = 'results/{dataset}/{population}/GatherEditing/AnnotatedSites.bed'
    output:
        SAMPLE_OUT = 'results/{dataset}/{population}/QTLprep/sample_name_list.txt',
        AP_raw_OUT = temp('results/{dataset}/{population}/QTLprep/scale_center/AP_raw.bed'),
        AP_norm_OUT = temp('results/{dataset}/{population}/QTLprep/ranknorm/AP_norm.bed'),
        EL_raw_OUT = temp('results/{dataset}/{population}/QTLprep/scale_center/EL_raw.bed'),
        EL_norm_OUT = temp('results/{dataset}/{population}/QTLprep/ranknorm/EL_norm.bed')
    params:
        MIN_N_SAMPLES = 10,
        MIN_DP_PER_SAMPLE = 10,
        MIN_AP_PER_SAMPLE = 1,
        MIN_AP_ROWSUM = 25,
        EL_ROWMEAN_RANGE = "0.001 0.99",
        EXCLUDE_SAMPLES = "NA18855" # this sample is an anomaly
    threads: 1
    resources: time=2000, mem_mb=25000, cpu=1
    script: '../scripts/prep_phenoBED_v2.R'


rule Index_phenotype:
    message: '### Index phenotype BED file'
    input: 
        EL_norm = rules.prep_phenotype.output.EL_norm_OUT,
        EL_raw = rules.prep_phenotype.output.EL_raw_OUT,
        AP_norm = rules.prep_phenotype.output.AP_norm_OUT,
        AP_raw = rules.prep_phenotype.output.AP_raw_OUT
    output: 
        EL_norm = 'results/{dataset}/{population}/QTLprep/ranknorm/EL_norm.bed.gz',
        EL_raw = 'results/{dataset}/{population}/QTLprep/scale_center/EL_raw.bed.gz',
        AP_norm = 'results/{dataset}/{population}/QTLprep/ranknorm/AP_norm.bed.gz',
        AP_raw = 'results/{dataset}/{population}/QTLprep/scale_center/AP_raw.bed.gz'
    shell:
        '''
        INPUTS=({input})
        OUTPUTS=({output})
        for i in {{0..3}}; do
            bgzip -c ${{INPUTS[i]}} > ${{OUTPUTS[i]}} 
            tabix -0 -s 1 -b 2 -e 3 ${{OUTPUTS[i]}}
        done
        '''
    

rule PermuteAndPCA_Phenotype:
    message: '### Run PCA on phenotype with permutation'
    input: 'results/{dataset}/{population}/QTLprep/ranknorm/EL_norm.bed.gz'
    output: 'results/{dataset}/{population}/QTLprep/pheno_permuted_pca.pca'
    shell:
        '''
        Rscript workflow/scripts/PermuteAndPCA.R {input} {output} 
        '''


rule prep_genotype:
    message: '### Prepare genotype vcf file'
    input: 
        sample_file = 'results/{dataset}/{population}/QTLprep/sample_name_list.txt',
        vcf =  '/project2/yangili1/zpmu/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz'
    output: 'results/{dataset}/{population}/Genotype/1kg_p3_genotype.{chrom}.vcf.gz'
    params:
        min_MAF = 0.05,
        max_HWE = 1e-3
    threads: 2
    resources: time=2000, mem_mb=15000, cpu=2
    shell:
        '''
        bcftools view \
            --threads {threads} \
            --samples-file {input.sample_file} \
            -i "AF >= {params.min_MAF} && HWE < {params.max_HWE}" \
            -Oz -o {output} \
            {input.vcf}
        
        bcftools index --threads {threads} --tbi {output}
        '''

rule ConcateGenotypeVCF:
    message: '### Concatenate Genotype VCF'
    input: expand('results/{{dataset}}/{{population}}/Genotype/1kg_p3_genotype.{chrom}.vcf.gz', chrom=AUTOSOMES)
    output: 'results/{dataset}/{population}/Genotype/1kg_p3_genotype_chr1-22.vcf.gz'
    resources: cpu = 1, mem = 25000, time = 1000
    shell: 
        '''
        bcftools concat {input} -O z -o {output}
        bcftools index {output}
        '''

rule pca_genotype:
    message: '### Run pca on genotype'
    input: 'results/{dataset}/{population}/Genotype/1kg_p3_genotype_chr1-22.vcf.gz'
    output: 'results/{dataset}/{population}/QTLprep/geno_pca.pca'
    params:
        out_prefix = 'results/{dataset}/{population}/QTLprep/geno_pca'
    shell:
        '''
            QTLtools pca \
            --seed 123 \
            --maf 0.01 \
            --vcf {input}  \
            --out {params.out_prefix} \
            --center \
            --scale
        '''


rule MakeCovarianceMatrix_FixedPCs:
    '''
    Here fixed PCs means I used fixed number of genotype pca PCs.
    For phenotype, the number of PCs used are more varied. Details see script.
    '''
    message: '### Make covariance matrix with fixed PCs'
    input:
        PhenoPCs = 'results/{dataset}/{population}/QTLprep/pheno_permuted_pca.pca',
        GenoPCs = 'results/{dataset}/{population}/QTLprep/geno_pca.pca'
    output: 'results/{dataset}/{population}/QTLprep/CovarianceMatrix_FixedPCs.txt'
    params:
        Fixed_Geno_PCs = 4
    shell:
        '''
            cat {input.PhenoPCs} <( awk 'NR>1' {input.GenoPCs} | head -{params.Fixed_Geno_PCs} ) 1> {output}
        '''




#----------------------------------------------------------------------------
# SECTION 2: MAP QTLS (PERMUTATION PASS)
#----------------------------------------------------------------------------


rule MapQTL_FixedPCs:
    input:
        vcf = 'results/{dataset}/{population}/Genotype/1kg_p3_genotype.{chrom}.vcf.gz',
        bed = 'results/{dataset}/{population}/QTLprep/ranknorm/EL_norm.bed.gz',
        cov = 'results/{dataset}/{population}/QTLprep/CovarianceMatrix_FixedPCs.txt'
    output: 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_{window}/{chrom}/permutations_pass.txt'
    params:
        cis_window = '{window}'
    resources: cpu = 1, mem = 12000, time = 1000
    shell:
        '''
        # already ranknorm normalized, do not the normal option 
        QTLtools cis \
            --seed 123 --silent \
            --vcf {input.vcf} --bed {input.bed} --cov {input.cov}  --out {output} \
            --window {params.cis_window} \
            --permute 1000 \
            --region {wildcards.chrom}
        '''


rule GatherQTLtoolsOutput_FixedPC_Perm:
    input: expand('results/{{dataset}}/{{population}}/edQTL_FixedPCs/ranknorm/cis_{{window}}/{CHROM}/permutations_pass.txt', CHROM=AUTOSOMES)
    output: 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_{window}/permutations_pass_all_chr.txt.gz'
    log: 'logs/GatherQTLtoolsOutput_FixedPC/{dataset}_{population}_cis_{window}_Perm.log'
    shell:
        '''
        (cat {input} | gzip - > {output}) &> {log}
        '''


rule AddQvalueToPermutationPass:
    '''
    Add qvalue
    '''
    message: '### Add qvalue to the output of QTLtools permutation pass'
    input: rules.GatherQTLtoolsOutput_FixedPC_Perm.output
    output: 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_{window}/permutations_pass_all_chr.AddedQvalue.txt.gz'
    log: 'logs/AddQvalueToPermutationPass/{dataset}_{population}_cis_{window}.log'
    shell:
        '''
            Rscript workflow/scripts/AddQvalueToQTLtoolsOutput.R {input} {output}
        '''




#----------------------------------------------------------------------------
# SECTION 3: MAP QTLS (PERMUTATION) USING RANDOMIZED PHENOTYPE AS CONTROL
#----------------------------------------------------------------------------


rule PermutePhenoBed:
    '''    
    Permute the columns of phenotypes
    '''
    input: 'results/{dataset}/{population}/QTLprep/ranknorm/EL_norm.bed.gz'
    output: 'results/{dataset}/{population}/Permuted/ranknorm/EL_norm_permuted.bed.gz'
    params:
        prefix = 'results/{dataset}/{population}/Permuted/ranknorm/EL_norm_permuted'
    shell:
        '''
        Rscript workflow/scripts/PermutePhenotypeColumns.R {input} {params.prefix}.bed
        bgzip {params.prefix}.bed
        tabix -p bed {output}
        '''

rule PermuteAndPCA_Phenotype_ShuffleSamples:
    input: rules.PermutePhenoBed.output
    output: 'results/{dataset}/{population}/Permuted/ranknorm/pheno_permuted_pca.pca'
    shell: 'Rscript workflow/scripts/PermuteAndPCA.R {input} {output}'


use rule MakeCovarianceMatrix_FixedPCs as MakeCovarianceMatrix_PermutedSamples with:
    input:
        PhenoPCs = rules.PermuteAndPCA_Phenotype_ShuffleSamples.output,
        GenoPCs = 'results/{dataset}/{population}/QTLprep/geno_pca.pca'
    output: 'results/{dataset}/{population}/Permuted/ranknorm/CovarianceMatrix_FixedPCs.txt'


use rule MapQTL_FixedPCs as mapQTL_Permuted with:
    input:
        vcf = 'results/{dataset}/{population}/Genotype/1kg_p3_genotype.{chrom}.vcf.gz',
        bed =  rules.PermutePhenoBed.output, # use permuted phenotype
        cov = rules.MakeCovarianceMatrix_PermutedSamples.output
    output: 'results/{dataset}/{population}/Permuted/ranknorm/cis_{window}/{chrom}/permutations_pass.txt'


rule GatherQTLtoolsOutput_PermutedSamples:
    input: expand('results/{{dataset}}/{{population}}/Permuted/ranknorm/cis_{{window}}/{CHROM}/permutations_pass.txt', CHROM=AUTOSOMES)
    output: 'results/{dataset}/{population}/Permuted/ranknorm/cis_{window}/permutations_pass_all_chr.txt.gz'
    shell:
        '''
        cat {input} | gzip - > {output}
        '''


rule AddQvalueToPermutationPass_PermutedSamples:
    '''
    Add qvalue
    '''
    message: '### Add qvalue to the output of QTLtools permutation pass'
    input: rules.GatherQTLtoolsOutput_PermutedSamples.output
    output: 'results/{dataset}/{population}/Permuted/ranknorm/cis_{window}/permutations_pass_all_chr.AddedQvalue.txt.gz'
    shell:
        '''
            Rscript workflow/Scripts/AddQvalueToQTLtoolsOutput.R {input} {output}
        '''




#----------------------------------------------------------------------------
# SECTION 4: MAP QTLs (nominal pass) 
#----------------------------------------------------------------------------

rule MapQTL_Nom_FixedPCs:
    """
    output columns:
    1. phenotype_id
    2. phenotype_chr
    3. phenotype_start
    4. phenotype_end
    5. phenotype_strand
    6. num_variant_tested
    7. variant_dist
    8. genotype_id
    9. genotype_chr
    10. genotype_start
    11. genotype_end
    12. pval_nominal
    13. pval_r2
    14. slope
    15. top_variant_flag
    """
    message: '### Run QTLtools nominal pass with fixed PCs'
    input: 
        vcf = 'results/{dataset}/{population}/Genotype/1kg_p3_genotype.{chrom}.vcf.gz',
        bed = 'results/{dataset}/{population}/QTLprep/ranknorm/EL_norm.bed.gz',
        cov = 'results/{dataset}/{population}/QTLprep/CovarianceMatrix_FixedPCs.txt'
    output: temp('results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_{window}/{chrom}/Nominal.txt')
    resources: cpu = 1, mem = 15000, time = 60
    shell: 
        '''
        QTLtools cis \
            --seed 123 --silent \
            --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --out {output} \
            --nominal 1 \
            --normal --window {wildcards.window} \
            --region {wildcards.chrom}
        '''

rule TabixNominalPassOutFile:
    message: '### bgzip and tabix nominal output'
    input: rules.MapQTL_Nom_FixedPCs.output
    output: 'results/{dataset}/{population}/edQTL_FixedPCs/ranknorm/cis_{window}/{chrom}/Nominal.txt.gz'
    log: 'logs/TabixNominalPassOutFile/{dataset}_{population}_cis_{window}/{chrom}.log'
    shell: 
        '''
        (awk 'BEGIN {{OFS="\t"}}; {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}}' {input} | \
            sort -k9 -k10 -k11 -V | \
            bgzip > {output}) &> {log}
        (tabix -s 9 -b 10 -e 11 {output}) &> {log}
        '''




#----------------------------------------------------------------------------
# SECTION 5: EXTRA
#----------------------------------------------------------------------------

# test on permuted sample
# not part of main pipeline, so never add to rule all.


use rule MapQTL_Nom_FixedPCs as MapQTL_Nom_PermutedSamples with:
    message: '### Nominal pass on permuted samples'
    input: 
        vcf = 'results/{dataset}/{population}/Genotype/1kg_p3_genotype.{chrom}.vcf.gz',
        bed = 'results/{dataset}/{population}/Permuted/ranknorm/EL_norm_permuted.bed.gz',
        cov = 'results/{dataset}/{population}/Permuted/ranknorm/CovarianceMatrix_FixedPCs.txt'
    output: 'results/{dataset}/{population}/Permuted/ranknorm/cis_{window}/{chrom}/Nominal.txt'


