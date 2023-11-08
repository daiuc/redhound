# /usr/bin/env python

'''Prepare inputs for sashimi plots

Two main ingredients:
    - averaged bigwig files by allele type
    - linkage files 
'''

__author__    = "Chao Dai"
__email__     = "chaodai@uchicago.edu"
__status__    = "Development"
__version__   =  "v0.0.1"



import os
import gzip
import numpy as np
import pysam as ps
import pandas as pd
import pyBigWig as pw


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


def refineSamples(SamplesByGT, PhenoSamples, bwfiles):
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
        gt[k]['bwfiles'] = [x for x in bwfiles if any([s in x for s in v])]
    return gt


def avgBigwig(bwfiles, grange):
    '''Average multiple bigwig files in a specific region
    
    bwfiles : list of bigwig files (path).
    grange  : tuple. Genomic range, BED like 0-based coordinates, 
              eg. ('chr1', 25101, 27101)
    
    return  : a dictionary of keys: 
        - header, for pyBigWig to write as header
        - values, for pyBigwig to addentries as values
    '''
    chrom, start, end = grange
    values = []
    bwo = {}
    for bw in bwfiles:
        if not os.path.isfile(bw): continue
        with pw.open(bw, 'rt') as b:
            # header = list(b.chroms().items())
            header = [(chrom, b.chroms()[chrom])]
            vals = b.values(chrom, start, end, numpy=True)
            vals = np.nan_to_num(vals)
            values.append(vals)
    
    if values != [] and header != []:
        avgValues = np.mean(values, axis=0)
        bwo = {'header': header, 'values': avgValues}
    return bwo


def strTofloat(frac):
    '''Convert string fraction to float
    frac : str. eg. '5/10'
    return : float
    '''
    a, b = frac.split('/')
    a, b = int(a), int(b)
    if b == 0: frac = 0
    else: frac = a/b
    return frac 


def readPsiFromTabix(tabixfile, grange, pid):
    '''Read splicing PSI file (tabixed) based on region and pid
       tabixfile : str, bgziped and tabix indexed psi file
       grange    : tuple, plot region, eg. ('chr1', 1000, 2000)
       pid       : str, pid, used to filter only the intron desired
    '''
    
    psi_cols = pd.read_csv(tabixfile, sep='\t', nrows=2)
    psi_cols = [c.replace('#', '') for c in psi_cols.columns]
    psi = ps.TabixFile(tabixfile)
    
    df = []
    c, s, e = grange
    for row in psi.fetch(c, s, e, parser=ps.asTuple()):
        df.append(row[0:])
    df = pd.DataFrame(df, columns=psi_cols)    
    gid = df[df.pid == pid]['gid']
    df = df[df.gid.isin(gid)]
    return(df)


def makeLinkTable(PsiDF, samples):
    '''from PSI dataframe for a given intron cluster get link table for sashimi
       PsiDF : dataframe, PSI table limited to specific region and intron cluster
       samples: list, sample name list for a given genotype
    '''
    mean = PsiDF[samples].apply(pd.to_numeric, errors='coerce', axis=1).mean(axis=1)
    df = pd.DataFrame({'Chr1': PsiDF.Chr, 'start1': PsiDF.start, 'start2': PsiDF.start, 
                       'Chr2': PsiDF.Chr, 'end1': PsiDF.end, 'end2': PsiDF.end,
                       'Mean': mean})
    df = df[df['Mean'] > 1e-5]
    return df


def readA2IFromTabix(tabixfile, grange, pids, dropValues=True):
    '''Read A2I editing sites from a given range
    tabixfile : str, bgzipped and tabix indexed gathered Editing level file
    grange    : tuple, plot region, eg. ('chr1', 1000, 2000), note start and end
                are integers
    pids      : list, of A2I pids that are edQTL hits
    dropValues: bool, drop editing level value columns if true
    '''
    a2i_cols = pd.read_csv(tabixfile, sep='\t', nrows=2)
    a2i = ps.TabixFile(tabixfile)
    a2i_cols = a2i.header[0].split('\t')
    print(len(a2i_cols))
    df = []
    c, s, e = grange
    for row in a2i.fetch(c, s, e, parser=ps.asTuple()):
        df.append(row[0:])
    df = pd.DataFrame(df, columns=a2i_cols)
    if 'name' in a2i_cols:
        pidcol = 'name'
    elif 'pid' in a2i_cols:
        pidcol = 'pid'
    df = df[df[pidcol].isin(pids)]
    df.reset_index(inplace=True, drop=True)
    
    if dropValues:
        df = df.iloc[:, 0:6]
    return(df)



def writeIni(template, outfile, fdict, alleles, gtf, a2i):
    '''Generate ini file for plotting sashimi plots
    template : str. Path to a template ini file.
    outfile  : str. Path of ini file to write into.
    fdict    : dict. Dictionary hosting bw file path, and link file path
                eg {0: {'bw': 'allele0.bw', 'link': 'lnk0.tsv'}, 1:...}
    alleles  : dict. eg {0: '0_AA', 1: '1_AC', 2: '2_CC'}
    title    : str, title name,
    gtf      : str, gtf file path
    a2i      : bed file for a2i hits in the region
    
    return : None.
    '''
    
    from jinja2 import Template
    
    print(f'Generate ini files: {outfile} ...\n')
    with open(template) as t:
        t_content = t.read()
        tt = Template(t_content)
        outini = tt.render( # out ini
                # title = title,
                # title_height = 1,
                where = "right",
                title_fontsize = 4,
                allele0 = alleles[0] if len(alleles) == 3 else '0',
                bw0 = fdict[0]['bw'],
                lnk0 = fdict[0]['link'],
                allele1 = alleles[1] if len(alleles) == 3 else '1',
                bw1 = fdict[1]['bw'],
                lnk1 = fdict[1]['link'],
                allele2 = alleles[2] if len(alleles) == 3 else '2',
                bw2 = fdict[2]['bw'],
                lnk2 = fdict[2]['link'],
                gtf_file = gtf,
                edQTL_bed = a2i
        )
    with open(outfile, 'w') as f:
        f.write(outini + '\n')

def main(options):
    
    outdir, outbase = os.path.dirname(options.outprefix), os.path.basename(options.outprefix)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    # bigwig files
    bwfiles = []
    with open(options.bwfiles) as f:
        bwfiles = [x.strip() for x in f.readlines()]

    # indexed bcf file for extracting genotypes
    vcf = ps.VariantFile(options.vcf)

    # here pid and snp are both for splicing
    pid = options.pid
    # standardize SNP and construct allele types
    snpid = options.snp
    snp = fixSNPformat(options.snp)
    if len(snp) == 4:
        alleles = {
            0: f'0_{snp[2]*2}',
            1: f'1_{snp[2] + snp[3]}',
            2: f'2_{snp[3]*2}'}
    else:
        alleles = {}
    
    # a2I editing pids with beta_adj pval < 0.01
    with open(options.a2iPids) as f:
        a2ipids = [x.strip() for x in f.readlines()]
    
    # plot range, use extend to extend both ends
    a, b = options.region.split(':')
    b, c = b.split('-')
    extend = int(options.extend)
    grange = a, int(b), int(c)
    
    # sample list in phenotype
    with open(options.samps) as f:
        pSamps = [l.strip() for l in f.readlines() if l.strip() != '']

    # get sample names and bw files
    samples = getSamplesByAllele(vcf, snpid)
    samples = refineSamples(samples, pSamps, bwfiles)
    
    #------------ operate on genotype 0, 1, and 2 ---------------------

    # PSI dataframe
    psidf = readPsiFromTabix(options.psiFile, grange, pid)

    avgBW = {} # avg bigwig for an allele
    outDict = {} # output file for each allele 0, 1, 2
    for k, d in samples.items():
        chrom, start, end = grange
        
        # ------------------ average bigwig ------------------------
        avgBW = avgBigwig(d['bwfiles'], (chrom, start-extend, end+extend))
        bwOut = os.path.join(outdir, f'{outbase}_{str(k)}.bw')

        print(f'Extracting average bigwig for allele {k}, write to {bwOut}...\n')
        with pw.open(bwOut, 'wt') as bw:
            bw.addHeader(avgBW['header'])
            bw.addEntries(chrom, start - extend, values=avgBW['values'], span=1, step=1)


        # ----------------- prep link table for splicing ------------

        # Extract PSI for allele and write link table, PSI is tabixed
        psi_cols = pd.read_csv(options.psiFile, sep='\t', nrows=2)
        psi_cols = psi_cols = [c.replace('#', '') for c in psi_cols.columns]

        lnkOut = os.path.join(outdir, f'{outbase}_links_{str(k)}.tsv')
        print(f'Extract link tables for allele {k}, write to {lnkOut}...\n')

        lnkdf = makeLinkTable(psidf, d['samples'])
        lnkdf.to_csv(lnkOut, sep='\t', header=False, index=False)
                
        # store written bw and link file paths by allele for ini generation, note only basename for bw file in ini
        outDict[k] = {'bw': os.path.basename(bwOut), 
                      'link': os.path.basename(lnkOut)}
        
        iniOut = os.path.join(outdir, f'{outbase}.ini')


        # ----------- prep A2I editing level by genotype ------------
        # elOut = os.path.join(outdir, f'{outbase}_EL_{str(k)}.tsv')
        # a2i_full = readA2IFromTabix(options.a2iFile, grange, a2ipids,
        #                          dropValues=False)
        # samplesIN = list(set(d['samples']).intersection(list(a2i_full.columns)))
        # eldf = pd.concat([a2i_full.iloc[:, :6], a2i_full.loc[:, samplesIN]],
        #                    axis=1)
        # eldf.to_csv(elOut, sep='\t', header=True, index=False)
    

    # ----------------- prep A2I site BED file ------------
    # prepare BED file for edQTL hits (editing sites) in this window
    print('Extract editing sites BED file for edQTL with beta-pval < 0.01.\n')
    a2i = readA2IFromTabix(options.a2iFile, grange, a2ipids)
    a2iOut = os.path.join(outdir, f'{outbase}_a2i.bed')
    a2i.to_csv(a2iOut, sep='\t', header=False, index=False)



    # ---------------- Generate ini file for sashimi plot -------------
    if options.template:
        writeIni(template=options.template, outfile=iniOut, 
                 fdict=outDict, alleles=alleles, gtf=options.gtf, a2i=a2iOut)

    print('Done.')


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-V", "--vcf", dest="vcf",
        type=str, required=True,
        help="gzipped and tbi indexed vcf file with sample genotypes")

    parser.add_argument("-B", "--bwfiles", dest="bwfiles",
        type=str, required=True,
        help="a file that includes a list of Bigwig file paths")
    
    parser.add_argument("--psiFile", dest="psiFile",
        type=str, required=True, help="phenotype of psi file")

    parser.add_argument("--a2iFile", dest="a2iFile",
        type=str, required=True, help="tabixed phenotype of gathered editing level file")

    parser.add_argument("-S", "--snp", dest="snp",
        type=str, required=True,
        help="SNP ID, as in the ID field in vcf, eg. chr1:1000:A:C, 1-based.")
    
    parser.add_argument("--pid", dest="pid",
        type=str, required=True,
        help="pid of intron, eg. 1:1489274:1490257:clu_2536_+")
    
    parser.add_argument("--a2iPids", dest="a2iPids",
        type=str, required=True,
        help="a file storing list of A2I pids that have beta adj. pval<0.01")
    
    parser.add_argument("-P", "--samps", dest="samps",
        type=str, required=True,
        help="A text file with 1 sample ID per line, indicating samples that \
              exist in phenotypes.")
    
    parser.add_argument("-R", "--region", dest="region",
        type=str, required=True,
        help="Phenotype genomic region, format: chr1:1000-2000, 0-base as in bed.")
    
    parser.add_argument("--gtf", dest="gtf",
        type=str, required=True,
        help="GTF file, e.g. Gencode primary annotation gtf.")
    
    parser.add_argument("-E", "--extend", dest="extend",
        default = 50, 
        help="Extend region on both ends to help plotting.")
    
    parser.add_argument("-T", "--template", dest="template",
        type=str, required=False, help="ini template")

    parser.add_argument("-O", "--outprefix", dest="outprefix",
        type=str, required=True, help="Output prefix, eg. path/to/file_prefix")
    
    options = parser.parse_args()
    
    main(options)