import pandas as pd
import numpy as np
from itertools import chain

# Hi-C utilities imports:
import cooler
import bioframe
import cooltools
import os
import argparse

parser = argparse.ArgumentParser(description='Call loops using a round donut, with cluster filtering.')
parser.add_argument('-i',required=True, metavar='mcool_dir',
                    help='.mcool file')
parser.add_argument('-e', required=True, metavar='expected',
                    help='expected file')
parser.add_argument('-r', required=True, metavar='resolution',
                    help='resolution for loop calling')

args = parser.parse_args()

binsize = int(args.r)

res = int(binsize/1000)
out = "%s_%sk_loops_rounddonut" % (os.path.basename(args.i).split('.')[0],str(res))


# Open cool file with Micro-C data:
clr = cooler.Cooler(f'./{args.i}::/resolutions/{binsize}')
expected = pd.read_table(args.e)


# # Use bioframe to fetch the genomic features from the UCSC.
# mm9_chromsizes = bioframe.fetch_chromsizes('mm9')
# mm9_cens = bioframe.fetch_centromeres('mm9')
# # create a view with chromosome arms using chromosome sizes and definition of centromeres
# mm9_arms = bioframe.make_chromarms(mm9_chromsizes, mm9_cens)

# # get rid of chromarm in name. bc mouse chromosomes acrocentric
# mm9_arms['name'] = mm9_arms['name'].apply(lambda x: x.rstrip('_p'))
# mm9_arms = mm9_arms[mm9_arms['name']!='chrM']
# mm9_arms = mm9_arms[mm9_arms['name']!='chrY']


# manually create list because bioframe is glitchy...
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
ends = [197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172, 129993255, 121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430, 166650296]

mm9_arms = pd.DataFrame()
mm9_arms['chrom'] = chroms
mm9_arms['start'] = 0
mm9_arms['end'] = ends 
mm9_arms['name'] = chroms


kernels = cooltools.api.dotfinder.recommend_kernels(10_000)
# create a grid of coordinates from -5 to 5, to define round kernels
# see https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html for details
half = 5  # half width of the kernel
x, y = np.meshgrid(
    np.linspace(-half, half, 2*half + 1),
    np.linspace(-half, half, 2*half + 1),
)
# now define a donut-like mask as pixels between 2 radii: sqrt(7) and sqrt(30):
mask = (x**2+y**2 > 7) & (x**2+y**2 <= 30)
mask[:,half] = 0
mask[half,:] = 0

# lowleft mask - zero out neccessary parts
mask_ll = mask.copy()
mask_ll[:,:half] = 0
mask_ll[half:,:] = 0

# new kernels with more round donut and lowleft masks:
kernels_round = {'donut': mask,
 'vertical': kernels["vertical"].copy(),
 'horizontal': kernels["horizontal"].copy(),
 'lowleft': mask_ll}

 #### call dots using redefined kernels (without clustering)

dots_round_df_all = cooltools.dots(
    clr,
    expected=expected,
    view_df=mm9_arms,
    kernels=kernels_round, # provide custom kernels
    max_loci_separation=2_000_000,
    clustering_radius=20_000,
    nproc=16,
    lambda_bin_fdr=0.05,
    n_lambda_bins=50
)


cols = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']

# output loops as a bedpe file and as a csv
dots_round_df_all.to_csv(out + '.csv', index=None)
dots_round_df_all[cols].to_csv(out + '.bedpe', index=None, header=None, sep='\t')





