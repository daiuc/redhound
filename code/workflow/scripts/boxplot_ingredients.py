#!/usr/bin/env python

'''Prepare input data for plotting genotype specific box plots for molecular 
QTLs.

Needed data:

- eQTL-edQTL coloc summary file: to get lead SNP
- VCF file to find sample names for Genotype 0, 1, 2.
- eQTL expression data for the colocalized phenotype for 3 genotypes
- edQTL editing level for the colocalized phenotype for 3 genotypes
'''

__author__    = "Chao Dai"
__email__     = "chaodai@uchicago.edu"
__status__    = "Development"
__version__   =  "v0.0.1"


import sys
import os
import gzip
import re
import json

import pandas as pd
import pysam as ps


#-------------------------- CustomFunctions -----------------------------------

def fixSNPformat(snp):
    '''conform SNP format to a tuple, e.g. ('chr1', 1000)
    snp : str, concatenated by ":"
    '''
    snp = snp.strip().split(':')
    if 'chr' not in snp[0]:
        snp[0] = 'chr' + snp[0]
    return tuple(snp)

def getIndex(l, m):
    '''returns index of all occurances of v in l.
    l : list to query from
    m : value match
    '''
    return [k for k,v in enumerate(l) if v == m]


def getSamplesByAllele(bcf, snpid):
    '''Get sample names by allele
    bcf   : pysam VariantFile object
    snpid : snp_id, e.g: 1:1000:G:A
    
    return : dict
        key   : 0, 1, or 2, representing allele
        value : sample names that belong to each allele
    '''
    snp = snpid.split(':')
    chrom, pos = str(snp[0]), int(snp[1])
    if not chrom.startswith('chr'):
        chrom = 'chr'+chrom
    
    samples = list(bcf.header.samples) # sample IDs
    records = []
    for rec in bcf.fetch(chrom, pos - 1, pos + 1):
        if rec.id == snpid:
            records.append(rec)
        elif rec.pos == pos: # in case snpid not matching
            records.append(rec)
    try:
        rec = records[0] # tolerate multipe records per pos or snpID
        gt = [x['GT'] for x in rec.samples.values()]
        gt = [a+b for a,b in gt]
        gt0 = getIndex(gt, 0)
        gt1 = getIndex(gt, 1)
        gt2 = getIndex(gt, 2)
    except IndexError:
        print(f'Error! No record found for this SNP!\n')
        exit(1)
    try:
        samples_by_allele = {
            0:[samples[x] for x in gt0],
            1:[samples[x] for x in gt1],
            2:[samples[x] for x in gt2]}
    except UnboundLocalError:
        print(f'Error! SNP is not found in vcf file!\n')
        exit(1)
    
    return samples_by_allele


def refineSamples(SamplesByGT, PhenoSamples):
    '''Ensure samples appear in both genotype and phenotypes
    --------------
    SamplesByGT  : dict. Dictionary of genotype and sample names, as in output
                   of getSamplesByAllele.
    PhenoSamples : list. List of samples that exist in phenotype, for instance
                   splicing. Sample names in genotype dict must be also in this
                   list.
    bwfiles      : list. File path of bigwigs.
    --------------
    return : dictionary of {0: {'sample':[], 'bwfiles':[]}}, where 0,1,2 is GT.
    '''
    gt = {}
    for k, v in SamplesByGT.items(): # k=0,1,2; v=[sampleIDs]
        gt[k] = {}
        v = [x for x in v if x in PhenoSamples]
        gt[k]['samples'] = v
    return gt

def readPhenoBed(tabixfile, grange, pid):
    '''Read editing Level PSI file (tabixed) based on region and pid
       tabixfile : str, bgziped and tabix indexed, qtl pheno bed file format.
                   i.e. first 6 cols like BED, 7+ are samples
       grange    : tuple, plot region, eg. ('chr1', 1000, 2000)
       pid       : str, pid, used to filter only the intron desired
    '''
    pheno = ps.TabixFile(tabixfile)
    header = pheno.header[0].split('\t')
    header[0:5] = ['#chr', 'start', 'end', 'pid', 'gid']
    
    df = []
    c, s, e = grange
    for row in pheno.fetch(c, s, e, parser=ps.asTuple()):
        df.append(row[0:])
    df = pd.DataFrame(df, columns=header)
    df = df[df.pid.isin([pid])]
    return(df)


def prepPhenoData(snpid, pid, phenoWindow, phenoFile, vcfFile):
    '''provide a snpid and phenotype id, get data for plotting the phenotype
       by genotype (0,1,2)
    '''
    snp = fixSNPformat(snpid)
    alleles = {0: snp[2]*2, 1: snp[2] + snp[3], 2: snp[3]*2}
    c, s, e = re.split('[:-]', phenoWindow)[:3]
    s, e = int(s) - 1, int(e)
    phenoSamples = ps.TabixFile(phenoFile).header[0].split('\t')[6:]
    bcf = re.sub('_chr\d\.', f'_{c}.', vcfFile)
    bcf = ps.VariantFile(bcf)
    
    # make sure pids from coloc summary match up with phenotype bed files
    if 'clu_' in pid: # for splicing pid
        pid = pid.split(':')[:4]
        pid = ':'.join(pid)
    elif 'ENSG' or 'ncRNA' in pid: # for eQTL pid or ncRNA pid
        pid = pid.split(':')[0]

    sampleByGenotype = getSamplesByAllele(bcf, snpid)
    sampleByGenotype = refineSamples(sampleByGenotype, phenoSamples)

    pheno = readPhenoBed(phenoFile, (c,s,e), pid)

    phenoByGenotype = {}
    for gt, d in sampleByGenotype.items():
        keepCols = ['#chr', 'start', 'end', 'pid', 'gid', 'strand'] + d['samples']
        phenoByGenotype[gt] = pheno[keepCols].melt(
            id_vars=keepCols[:6], var_name='Sample', value_name='measure'
        ).to_dict('list')

    jsonOut = {'pid': pid, 'snpid': snpid, 'GT': alleles, 'data': phenoByGenotype}
    return jsonOut






def main(options):

    from datetime import datetime

    colocFile = options.colocFile
    colocIdsFile = options.colocIdsFile
    phenoFile = options.phenoFile
    phenoName = options.phenoName # which phenotype to extract
    vcfFile = options.vcfFile
    outFile = options.outFile

    print(f'Started script at: {datetime.now()}')
    print(f'Extract hit QTLs based on coloc summary: {colocFile}, {colocIdsFile}')
    print(f'Extract phenotype data from {phenoFile}')
    print(f'Extract genotypes from {vcfFile}, with appropriate chromsome.')

    coloc = pd.read_csv(colocFile, sep='\t')
    coloc['id'] = [f'{x.strip()}|{y.strip()}' \
                   for x,y in zip(coloc.pheno1, coloc.pheno2)]

    with open(colocIdsFile) as f:
        colocIds = [x.replace(' ', '').strip() for x in f.readlines()]
    
    if all([x in list(coloc['id']) for x in colocIds]):
        coloc = coloc[coloc['id'].isin(colocIds)]
        coloc.sort_values('PP.H4.abf', ascending=False, inplace=True)
        coloc.reset_index(drop=True, inplace=True)
    else:
        print(f'IDs do not map completely! Check {colocFile} and {colocIdsFile}.')
        exit(1)

    if 'EDIT' in phenoName.upper():
        print('Extract editing by genotype.')
        coloc = coloc[['pheno1', 'p1window', 'p1topsnp']].drop_duplicates()
    else:
        print('Extract splicing/expression/ncRNA-expression by genotype.')
        coloc = coloc[['pheno2', 'p2window', 'p2topsnp']].drop_duplicates()
    coloc.columns = ['pheno', 'window', 'topsnp']

    jsonOut = {}
    for snpid, pid, window in zip(coloc['topsnp'], coloc['pheno'], coloc['window']):
            k = f'{pid}|{snpid}'
            jsonOut[k] = prepPhenoData(snpid, pid, window, phenoFile, vcfFile)
    
    print(f'\n\nTotal phenotypes: {len(jsonOut.keys())}')
    print(f'Writing extracted phenotype by genotype data in json format to {outFile}')
    with open(outFile, 'w') as f:
        json.dump(jsonOut, f)

    print(f'Done. {datetime.now()}.')





#-------------------------- Procedures ------------------------------


if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser()

    # colocFile = options.colocFile
    # colocIdsFile = options.colocIdsFile
    # phenoFile = options.phenoFile
    # phenoName = options.phenoName # which phenotype to extract
    # vcfFile = options.vcfFile
    # outFile = options.outFile
    
    parser.add_argument("--colocFile", dest="colocFile",
        type=str, required=True, help="Coloc summary file.")

    parser.add_argument("--colocIdsFile", dest="colocIdsFile",
        type=str, required=True,
        help="Hand picked coloc Ids, used for filtering coloc results.")

    parser.add_argument("--phenoFile", dest="phenoFile",
        type=str, required=True,
        help="Tabix indexed phenotype bed file, like editing, splicing, etc.")
    
    parser.add_argument("--phenoName", dest="phenoName",
        type=str, required=True,
        help="Phenotype name, options: editing, splicing, expression, ncRNA")
    
    parser.add_argument("--vcfFile", dest="vcfFile",
        type=str,
        default='/project2/yangili1/zpmu/1kg_b38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr2.filtered.shapeit2-duohmm-phased.vcf.gz',
        help="VCF file to getting sample names by genotype. Fine to leave default.")

    parser.add_argument("-O", "--outFile", dest="outFile",
        type=str, required=True, help="json data to write out to.")

    options = parser.parse_args()
    
    main(options)












