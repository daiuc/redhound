### Estimate number of A>I editing sites per kb in human genome
### Using different annotations

import pandas as pd

# Glogals
GENCODE_GTF = '/project2/yangili1/cdai/genome_index/hs38/gencode.v38.primary_assembly.annotation.gtf.gz'



# Use pygtftk to get CDS from gencode gtf file
rule Extract_Gencode:
  # input: 'test.gtf'
  input: GENCODE_GTF
  output: 'results/misc/Edit_site_per_kb/annotations/{feature}.bed'
  wildcard_constraints:
      feature = 'CDS|UTR'
  conda: 'gtftk'
  shell:
      '''
      # use pygtftk to get CDS from gencode gtf file and output bed format file
      gtftk select_by_key -i {input} -k gene_type -v protein_coding | \
        gtftk select_by_key -k feature -v {wildcards.feature} \
            -b -m gene_id,gene_name | \
        sort -V -k1,1 -k2,2n | uniq > {output}
      '''



# Because some features do overlap (within gene), we need to merge them
rule Merge_features:
  input: rules.Extract_Gencode.output
  output: 'results/misc/Edit_site_per_kb/annotations/{feature}.merged.bed'
  wildcard_constraints:
      feature = 'CDS|UTR'
  run:
      from pybedtools import BedTool

      # Step 1: Read the BED file into a DataFrame
      bed_file = input[0]
      print(f'Input file: {bed_file}')

      column_names = ['chrom', 'start', 'end', 'gid', 'score', 'strand']
      df = pd.read_csv(bed_file, sep='\t', header=None, names=column_names)

      # Step 2: Group the DataFrame by 'gid'
      print('Grouping by gid and strand')
      grouped = df.groupby(['gid', 'strand'])

      # Step 3: Merge overlapping features within each group
      print('Merging features')
      merged_features = []
      for group_name, group_data in grouped:
          gid, strand = group_name
          sorted_group = group_data.sort_values(by=['start', 'end'])
          bedtool_group = BedTool.from_dataframe(sorted_group)
          merged_group = bedtool_group.merge().to_dataframe()
          merged_group['gid'] = gid
          merged_group['score'] = "."
          merged_group['strand'] = strand
          
          merged_features.append(merged_group)

      # Step 4: Combine the merged features into a single BED file
      print('Combining features')
      merged_bed = pd.concat(merged_features)

      # Write to output - tab delimited, no index, no header
      print(f'Output file: {output[0]}')
      merged_bed.to_csv(output[0], sep='\t', index=False, header=False)



# merge introns
rule Merge_introns:
  input: '/project2/yangili1/cdai/annotations/hg38/use_ucsc_table_browser/gencode.v43_intron_proteinOnly.bed'
  output: 'results/misc/Edit_site_per_kb/annotations/intron.merged.bed'
  run:
      from pybedtools import BedTool

      # Step 1: Read the BED file into a DataFrame
      bed_file = input[0]
      print(f'Input file: {bed_file}')

      column_names = ['chrom', 'start', 'end', 'gid', 'score', 'strand']
      df = pd.read_csv(bed_file, sep=',')
      df.columns = column_names

      # clean up gid and remove duplicates
      df['gid'] = df['gid'].str.split('|').str[0]
      df = df.drop_duplicates()

      # Step 2: Group the DataFrame by 'gid'
      print('Grouping by gid and strand')
      grouped = df.groupby(['chrom', 'gid', 'strand'])

      # Step 3: Merge overlapping features within each group
      print('Merging features')
      merged_features = []
      for group_name, group_data in grouped:
        chrom, gid, strand = group_name
        sorted_group = group_data.sort_values(by=['chrom', 'start', 'end'])
        bedtool_group = BedTool.from_dataframe(sorted_group)
        merged_group = bedtool_group.merge().to_dataframe()
        merged_group['gid'] = gid
        merged_group['score'] = "."
        merged_group['strand'] = strand

        merged_features.append(merged_group)

      # Step 4: Combine the merged features into a single BED file
      print('Combining features')
      merged_bed = pd.concat(merged_features)

      # Write to output - tab delimited, no index, no header
      print(f'Output file: {output[0]}')
      merged_bed.to_csv(output[0], sep='\t', index=False, header=False)




# infer 5' UTR and 3' UTR from CDS and UTR annotations from gencode
# 5' UTR are UTRs upstream of CDS if on + strand, or downstream of CDS if on - strand
rule Infer_UTR_5p3p:
  input: 
    utr = 'results/misc/Edit_site_per_kb/annotations/UTR.merged.bed',
    cds = 'results/misc/Edit_site_per_kb/annotations/CDS.merged.bed'
  output: 
    flag = touch('results/misc/Edit_site_per_kb/annotations/inferUTR.done'),
    utr5p = 'results/misc/Edit_site_per_kb/annotations/5pUTR.merged.bed',
    utr3p = 'results/misc/Edit_site_per_kb/annotations/3pUTR.merged.bed',
  params:
    out_prefix = 'results/misc/Edit_site_per_kb/annotations',
    R_script = 'workflow/analysis-scripts/infer_5p_3p_utr.R'
  threads: 6
  shell:
    '''
    Rscript {params.R_script} {input.utr} {input.cds} {params.out_prefix} {threads}
    '''


# sort features by chromosome and start position
# also restrain chromosomes to 1-22, X, Y
rule Sort_features:
  input: 'results/misc/Edit_site_per_kb/annotations/{feature}.merged.bed'
  output: 'results/misc/Edit_site_per_kb/annotations/{feature}.merged.sorted.bed'
  wildcard_constraints:
      feature = 'CDS|UTR|intron|[53]pUTR'
  shell:
      '''
      sort -V -k1,1 -k2,2n {input} | \
        awk 'BEGIN {{OFS="\t"}}; $1 ~ /^chr[0-9]+$|^chrX$|^chrY$/ {{print $0}}' > {output}
      '''


# Normalize coordinates (still like BED format 0-based) such that within each
# (chr, strand) group, coordinates start from 0, and are tiled continuously.
# e.g. [0, 50), [50, 100), [100, 150), ...]
rule Normalize_coordinates:
  input: rules.Sort_features.output
  output: 'results/misc/Edit_site_per_kb/annotations/{feature}.merged.sorted.normed.tsv'
  wildcard_constraints:
      feature = 'CDS|UTR|intron|[53]pUTR',
      dataset = 'mRNA|chRNA'
  params:
    py_script = 'workflow/analysis-scripts/normalize_coordinates.py',
  shell:
      '''
      python {params.py_script} {input} {output}
      '''


# compute number of editing sites per window
rule Count_editing_per_window:
  input: 
    anno = rules.Normalize_coordinates.output,
    edit = 'results/{dataset}/YRI/GatherEditing/EL.txt'
  output: 
    flag = touch('results/misc/Edit_site_per_kb/count_edits_per_window/{feature}_{dataset}_YRI_countSitesPerWindow.done'),
    edit = 'results/misc/Edit_site_per_kb/count_edits_per_window/{feature}_{dataset}_YRI.editing-sites.norm.bed',
    counts = 'results/misc/Edit_site_per_kb/count_edits_per_window/{feature}_{dataset}_YRI.editing-per-window.count.bed'
  wildcard_constraints:
      feature = 'CDS|UTR|intron|[53]pUTR',
      dataset = 'mRNA|chRNA'
  params:
    wsize = 1000,
    out_prefix = 'results/misc/Edit_site_per_kb/count_edits_per_window/{feature}_{dataset}_YRI',
    R_script = 'workflow/analysis-scripts/count-editing-sites-per-window.R'
  shell:
    '''
    Rscript {params.R_script} -A {input.anno} -I {input.edit} -w {params.wsize} --out {params.out_prefix}
    ls {output.edit} {output.counts}

    '''


# # Use bedtools to get the intersection of editing sites and annotations.
# # report number of unique editing sites per annotation feature
# rule Count_editing_per_feature:
#   input: 
#     feature = rules.Sort_features.output,
#     editing = 'results/{dataset}/YRI/GatherEditing/EL.txt'
#   output: 'results/misc/Edit_site_per_kb/count_edits_to_features/{feature}_{dataset}_YRI.tsv'
#   wildcard_constraints:
#       feature = 'CDS|UTR|intron|infer.[53]p_utr',
#       dataset = 'mRNA|chRNA'
#   shell:
#       '''
#       bedtools intersect -a {input.feature} \
#         -b <(awk 'BEGIN {{OFS="\t"}}; NR > 1 {{print $1,$2,$3,$4,$5,$6}}' {input.editing}) \
#         -wa -c -s > {output}
#       '''





































