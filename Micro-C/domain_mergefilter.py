#!/usr/bin/env python
# coding: utf-8


# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Import python package for working with cooler files and tools for analysis
import cooler
import cooltools.lib.plotting
from cooltools import insulation

import glob 
import os.path

import pyBigWig

from pybedtools import BedTool 

import cooltools
cooltools.__version__

import bioframe
import argparse

parser = argparse.ArgumentParser(description='merge and filter domain calls for untreated/treated sample')
parser.add_argument('-m',required=True, metavar='mcool',
                    help='.mcool file')
parser.add_argument('-e', required=True, metavar='exp',
                    help='expected file')
parser.add_argument('-t', required=True, metavar='domain_dir',
                    help='expected file')
parser.add_argument('-n',required=True, metavar='mcool_treatment',
                    help='.mcool file')
parser.add_argument('-f', required=True, metavar='exp_treatment',
                    help='expected file')
parser.add_argument('-u', required=True, metavar='domain_dir_treatment',
                    help='expected file')
parser.add_argument('-o', required=True, metavar='domain_dir_treatment',
                    help='expected file')


args = parser.parse_args()

sample_dir = args.t
coolfile = args.m
expfile = args.e

sample_dirT = args.u
coolfileT = args.n
expfileT = args.f


outprefix = args.o

clr = cooler.Cooler(coolfile + '::resolutions/10000')
exp = pd.read_table(expfile)

clrT = cooler.Cooler(coolfileT + '::resolutions/10000')
expT = pd.read_table(expfileT)


# useful variables
chroms = list(range(1,20))
chroms.append("X")

bedcols = ['chrom', 'start', 'end']
bedcols2 = ['chr', 'start', 'end']




def merge_boundaries_in_tadlist(og):
    # flatten tad list to boundary list
    og2 = og[['chr', 'end']].copy()
    og2.columns = ['chr', 'start']

    og2 = og2.append(og[['chr', 'start']])

    # add an 'end'
    og2['end'] = og2['start'] + 1
    og2 = og2.sort_values(by='start').drop_duplicates()

    clr.chroms()[:].to_csv('chroms.txt', sep='\t', header=None, index=None)
    
    og2_bt = BedTool.from_dataframe(og2).sort(g='chroms.txt')

    merged = og2_bt.merge(d=80000)

    # intersect original boundaries with new ones
    cols = ['chr', 'oldstart', 'oldend', 'chrom2', 'merge1', 'merge2']
    new_bounds = og2_bt.intersect(merged, wa=True, wb=True).to_dataframe(names = cols)
    
    # calculate avg of merged boundaries
    new_bounds['merge2'] = new_bounds.merge2 - 1
    new_bounds['avg'] = new_bounds[['merge1', 'merge2']].apply(np.mean, axis=1)

    # redo domain list wiht new boudnaries
    new_bounds['start'] = new_bounds['oldstart']
    new_bounds = new_bounds[['chr', 'start', 'avg']]
    new_bounds.columns = ['chr', 'start', 'newstart']
    og = og.merge(new_bounds, on=['chr', 'start'], how='left')

    new_bounds.columns = ['chr', 'end', 'newend']
    og = og.merge(new_bounds, on=['chr', 'end'], how='left')
    og = og.drop_duplicates(subset=['chr', 'newstart', 'newend']) # this is the adjusted domain list
    og[['newstart', 'newend']] = og[['newstart', 'newend']].astype(int)

    return og


def generate_finetuned_domainlist(sample_dir, clr):
    # merge boundaries within 80kb
    final = pd.DataFrame() # final df for stroing tads for all chr

    for chrom in chroms:

        chrname = 'chr' + str(chrom)
        t1 = pd.read_csv(sample_dir + chrname + "_sweep1_hierTADS.csv")
        t1['chr'] = chrname
        t1 = t1[['chr', 'start', 'end']]

        t2 = pd.read_csv(sample_dir + chrname + "_sweep2_hierTADS.csv")
        t2['chr'] = chrname
        t2 = t2[['chr', 'start', 'end']]

        # read in tad lists for a chromosome
        og = t1.append(t2)

        og = merge_boundaries_in_tadlist(og)

        final = final.append(og)


    final['tadsize'] = final['newend'] - final['newstart']
    final = final[final['tadsize'] !=0]

    final = final[['chr', 'newstart', 'newend', 'tadsize']]
    
    # calculate insulation scores
    
    # calculate insulation score
    window = 120000
    df = insulation(clr, window, append_raw_scores=True,
                    verbose=True)

    bs_col = "boundary_strength_" + str(window)
    is_col = "log2_insulation_score_" + str(window)
    ib_col = "is_boundary_" + str(window)


    # finetune boundaries
    # flatten domain into boundaries again
    # flatten tad list to boundary list with +/- 60kb window
    final2 = final[['chr', 'newend']].copy()
    final2.columns = ['chr', 'newstart']
    final2 = final2.append(final[['chr', 'newstart']])
    final2['newstart'] = (final2['newstart'] - 60000).astype(int)

    # add an 'end'
    final2['newend'] = (final2['newstart'] + 120000).astype(int)

    # merge 
    insul_filtered = BedTool.from_dataframe(df[df.sum_counts_120000 > 12][bedcols + [is_col]].dropna())


    final2_bt = BedTool.from_dataframe(final2)

    cols = ['chrom', 'us', 'ds', 'chr', 'finestart', 'fineend', is_col]
    finetune = final2_bt.intersect(insul_filtered, wa=True, wb=True).to_dataframe(names = cols)

    # sort by insulation score - ascending- and keep the first/lowest log2 IS
    finetune = finetune.sort_values(by=is_col).drop_duplicates(subset=['chr','us','ds'])

    # take median of us and ds to get the original boundary
    finetune['merge_col'] = finetune[['us', 'ds']].apply(np.mean,axis=1).astype(int)
    
    # replace boundaries with fine tuned boundaries 
    #finestart was the start of the 10k insulation bin
    finetune = finetune[['chr', 'merge_col', 'finestart', 'log2_insulation_score_120000']]

    # merge us boundary
    final.columns = ['chr', 'merge_col', 'newend','tadsize']
    finaltads = final.merge(finetune, on=['chr', 'merge_col'], how='left')

    # merge ds boundary
    finaltads.columns = ['chr', 'drop1', 'merge_col','tadsize', 'finestart', 'log2_insulation_score_120000']
    finaltads = finaltads.merge(finetune, on=['chr', 'merge_col'], how='left', suffixes=("_us", '_ds'))
    
    
    # drop extra columns
    cols = ['chr','finestart_us', 'finestart_ds', is_col+'_us', is_col+'_ds']
    finaltads = finaltads[cols]

    # extra filtering locations that didnt have finetuned boundary/IS info
    finaltads = finaltads.dropna()

    # small spurious tads < 100kb
    finaltads = finaltads[(finaltads.finestart_ds - finaltads.finestart_us) >= 1e5]

    finaltads[['finestart_us', 'finestart_ds']] = finaltads[['finestart_us', 'finestart_ds']].astype(int)

    # get blacklist regions and remove domains overlapping with them
    # blacklist = <12 raw counts
    blacklist = BedTool.from_dataframe(df[df.sum_counts_120000 < 12][bedcols])

    finaltads_bt = BedTool.from_dataframe(finaltads)

    finaltads_filtered = finaltads_bt.intersect(blacklist, v=True).to_dataframe(names = cols)

    # merge outlier adjacent boundaries once more
    cols = ['chr', 'start', 'end', is_col + "_us", is_col+'_ds']
    finaltads_filtered.columns = cols
    finaltads_filtered = merge_boundaries_in_tadlist(finaltads_filtered)
    
    # only keep the fine tuned boundary locations
    finaltads_filtered = finaltads_filtered[['chr', 'newstart', 'newend']]
    cols = ['chr', 'start', 'end']
    finaltads_filtered.columns = cols

    return finaltads_filtered, df


# # get finetuned boundaries for -aux and +aux


tads0, insul0 = generate_finetuned_domainlist(sample_dir, clr)
tads4, insul4 = generate_finetuned_domainlist(sample_dirT, clrT)


# # merge auxin sample



# expand untreated boundary to +/- 80kb
tads0['chrom1'] = tads0['chr']
tads0['start_window1'] = tads0.start - 80_000
tads0['start_window2'] = tads0.start + 80_000

tads0['chrom2'] = tads0['chr']
tads0['end_window1'] = tads0.end - 80_000
tads0['end_window2'] = tads0.end + 80_000

cols = tads0.columns.to_list()
cols = cols[3:] + cols[:4]
tads0_bt = BedTool.from_dataframe(tads0[cols])

# convert auxin domain list into bedpe
# expand untreated boundary to +/- 80kba
tads4['chrom1'] = tads4['chr']
tads4['start_window1'] = tads4.start 
tads4['start_window2'] = tads4.start + 10_000

tads4['chrom2'] = tads4['chr']
tads4['end_window1'] = tads4.end 
tads4['end_window2'] = tads4.end + 10_000
cols = tads4.columns.to_list()
cols = cols[3:] + cols[:4]
tads4_bt = BedTool.from_dataframe(tads4[cols])

# keep unique auxin domains
tads4_unique = tads4_bt.pairtopair(tads0_bt, 
                                   **{'type': 'notboth','is': True})


# make list of boundaries in each 160kb window
tads0_bounds = tads0[cols[3:6] + cols[6:7] + cols[8:9]]
tads0_bounds.columns = cols[0:3] + cols[6:8]

tads0_bounds = tads0_bounds.append(tads0[cols[0:3] + cols[6:8]])
tads0_bounds_bt = BedTool.from_dataframe(tads0_bounds)

# intersect left anchor, LOJ = left outer join
tads4_unique = tads4_unique.intersect(tads0_bounds_bt, wa=True, wb=True, loj=True)



newcols = ['chr0', 'chr0b', 'wstart0', 'wend0', 'chr0c', 'finestart0']
tads4_unique = tads4_unique.to_dataframe(names = cols[:-1] + newcols)


tads4_unique.columns



# take control boundary. if it's null (-1) then take origianl boundary
tads4_unique['newstart'] = tads4_unique[['finestart0', 'start']].apply(np.max,axis=1)
tads4_unique = tads4_unique[['chrom2', 'end_window1','end_window2',
                            'chr', 'newstart', 'end']]
tads4_unique = BedTool.from_dataframe(tads4_unique)

# repeat for downstream boudnary
tads4_unique = tads4_unique.intersect(tads0_bounds_bt, wa=True, wb=True, loj=True)
tads4_unique = tads4_unique.to_dataframe(names=['chrom2', 'end_window1','end_window2',
                            'chr', 'newstart', 'end'] + newcols[1:])



tads4_unique['newend'] = tads4_unique[['finestart0', 'end']].apply(np.max,axis=1)
tads4_unique = tads4_unique[['chr','newstart','newend']].drop_duplicates()





tads4_unique.columns = ['chr', 'start', 'end']

domains = tads4_unique.append(tads0[['chr', 'start', 'end']])


# # filter based on inner/outer stripe counts



# convert domains list into bedpe for compatibility w/ pileup
df = domains.copy() 
df['chrom1'] = df.chr
df['chrom2'] = df.chr
df['start1'] = df.start - 5000
df['end1'] = df.start + 5000
df['start2'] = df.end - 5000
df['end2'] = df.end + 5000

# make masks defining the stripes from zhang 2021
i, j = (17,17)

inner_horizontal = np.zeros((35,35))
inner_vertical = np.zeros((35,35))
outer_horizontal = np.zeros((35,35))
outer_vertical = np.zeros((35,35))


inner_horizontal[i + 1, j-8:j-4]=1
inner_horizontal[i + 2, j-7:j-3]=1
inner_horizontal[i + 3, j-6:j-2]=1
inner_horizontal[i + 1, j-8:j-4]=1
inner_horizontal[i + 4, j-5:j-1]=1

inner_vertical[i + 1:i + 5, j-4] = 1
inner_vertical[i + 2:i + 6, j-3] = 1
inner_vertical[i + 3:i + 7, j-2] = 1 
inner_vertical[i + 4:i + 8, j-1] = 1

outer_horizontal[i-8, j-17:j-13] = 1
outer_horizontal[i-7, j-16:j-12] = 1
outer_horizontal[i-6, j-15:j-11] = 1
outer_horizontal[i-5, j-14:j-10]= 1

outer_vertical[i + 10:i + 14, j + 5]=1
outer_vertical[i + 11:i + 15, j + 6]=1
outer_vertical[i + 12:i + 16, j + 7]=1
outer_vertical[i + 13:i + 17, j + 8]=1

combined = (inner_horizontal+outer_horizontal+inner_vertical+outer_vertical) > 0

# get chromsizes for cooler pileup function
mm9_chromsizes = pd.read_table('chroms.txt', header=None)
mm9_chromsizes.columns = ['chrom', 'end']
mm9_chromsizes['start'] = 1
mm9_chromsizes = mm9_chromsizes[['chrom', 'start', 'end']]
mm9_chromsizes['name'] = mm9_chromsizes.chrom

mm9_regions = mm9_chromsizes[~mm9_chromsizes.chrom.isin(['chrY', 'chrM'])]

# Make pileup for each heatmap
# cooltools pileup function will determine center point of each anchor
# and then add/subract flank value to the center
# Make pileup for each heatmap
# cooltools pileup function will determine center point of each anchor
# and then add/subract flank value to the center
nzmax = 12
for klr in [clr, clrT]:
    stack = cooltools.pileup(klr, 
                    df[[ 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']], 
                    expected_df=exp, 
                    flank=170000,view_df=mm9_regions, nproc=16)

    stack = np.nan_to_num(stack)

    # check for pixels > 30 o/e
    mask = np.dstack([combined] * np.shape(stack)[2])
    stack1 = stack * mask
    df['no_outlierpixels'] = (np.nanmax(stack1, axis=(0,1)) < 30)

    # inner stripes must have > 10 nonzero pixels
    mask = np.dstack([inner_horizontal] * np.shape(stack)[2])
    stack1 = stack * mask
    df['innerh'] = np.count_nonzero(stack1, axis=(0,1)) > nzmax

    mask = np.dstack([inner_vertical] * np.shape(stack)[2])
    stack1 = stack * mask
    df['innerv'] = np.count_nonzero(stack1, axis=(0,1)) > nzmax

    # outer stripes must have > 10 nonzero pixels
    mask = np.dstack([outer_horizontal] * np.shape(stack)[2])
    stack1 = stack * mask
    df['outerh'] = np.count_nonzero(stack1, axis=(0,1)) > nzmax/2

    mask = np.dstack([outer_vertical] * np.shape(stack)[2])
    stack1 = stack * mask
    df['outerv'] = np.count_nonzero(stack1, axis=(0,1)) > nzmax/2

    df['keep'] = df[['no_outlierpixels', 'innerh', 'innerv','outerh', 'outerv']].apply(np.min,axis=1)

    df = df[df.keep]


# final filter of really small domains
df['tadsize'] = df.end - df.start 
df = df[df.tadsize > 100_000]



df=df[['chr','start','end']]


#################################### SAVE RESULTS  #################################### 

is_col = "log2_insulation_score_120000" 


# SAVE BIGWIG OF FILTERED INSULATION SCORES
out = os.path.basename(coolfile).split('.')[0]
bw = pyBigWig.open(out + "_log2insulationscore.bw", "w")
hd = list(zip(clr.chroms()[:].name, clr.chroms()[:].length))
bw.addHeader(hd)
bw.addEntries(insul0.chrom.to_list(), 
              insul0.start.to_list(), 
              ends=insul0.end.to_list(), 
              values=insul0[is_col].to_list())
bw.close()
insul0.to_csv(out + "_insulation_table.csv", index=None)

out = os.path.basename(coolfileT).split('.')[0]
bw = pyBigWig.open(out + "_log2insulationscore.bw", "w")
hd = list(zip(clrT.chroms()[:].name, clrT.chroms()[:].length))
bw.addHeader(hd)
bw.addEntries(insul4.chrom.to_list(), 
              insul4.start.to_list(), 
              ends=insul4.end.to_list(), 
              values=insul4[is_col].to_list())
bw.close()
insul4.to_csv(out + "_insulation_table.csv", index=None)


# ANNOTATE DOMAINS AS 0H, 4H OR both
shared = tads4_bt.pairtopair(tads0_bt, 
                                   **{'type': 'both','is': True})

shared = shared.to_dataframe(header=None)
shared = shared[shared.columns[-4:-1]]
shared.columns = ['chr', 'start', 'end']
shared['domaincall_map'] = 'both'

tads4_unique['domaincall_map'] = 'auxin'

df = df.merge(shared.append(tads4_unique), 
              how='left', on=['chr', 'start', 'end'])
df['domaincall_map'] = df.domaincall_map.fillna('untreated')

# SAVE CSV and BED FILE WITH DOMAINS
df = df.drop_duplicates()
df.to_csv(outprefix + '_filteredTADS.csv', index=None)

df['filler'] = '.' #coolbox requires 6 columns to view tads
df['filler2'] = '.'


df.to_csv(outprefix + '_filteredTADS.bedpe', index=None, header=None, sep='\t')

out = os.path.basename(coolfile).split('.')[0]
df[df.domaincall_map!='auxin'].to_csv(out + '_filteredTADS.bedpe', 
                                      index=None, header=None, sep='\t')

out = os.path.basename(coolfileT).split('.')[0]
df[df.domaincall_map!='untreated'].to_csv(out + '_filteredTADS.bedpe', 
                                      index=None, header=None, sep='\t')

# SAVE BED FILE WITH BOUNDARIES 
bounds = pd.DataFrame()
bounds['chr'] = df['chr']
bounds['start'] = df['end']
bounds = bounds.append(df[['chr', 'start']])
bounds['end'] = (bounds['start'] + 1e4).astype(int)
bounds.drop_duplicates().to_csv(outprefix + '_filteredBOUNDARIES.bed', index=None, header=None, sep='\t')

# save boundary file for each sample
df2 = df.copy()

out = os.path.basename(coolfile).split('.')[0]
df = df2[df2.domaincall_map != 'auxin']
bounds = pd.DataFrame()
bounds['chr'] = df['chr']
bounds['start'] = df['end']
bounds = bounds.append(df[['chr', 'start']])
bounds['end'] = (bounds['start'] + 1e4).astype(int)
bounds.drop_duplicates().to_csv(out + '_filteredBOUNDARIES.bed', index=None, header=None, sep='\t')


out = os.path.basename(coolfileT).split('.')[0]
df = df2[df2.domaincall_map != 'untreated']
bounds = pd.DataFrame()
bounds['chr'] = df['chr']
bounds['start'] = df['end']
bounds = bounds.append(df[['chr', 'start']])
bounds['end'] = (bounds['start'] + 1e4).astype(int)
bounds.drop_duplicates().to_csv(out + '_filteredBOUNDARIES.bed', index=None, header=None, sep='\t')

