# import core packages
import warnings
warnings.filterwarnings("ignore")
from itertools import combinations

# import semi-core packages
import matplotlib.pyplot as plt
from matplotlib import colors
plt.style.use('seaborn-poster')
import numpy as np
import pandas as pd

# import open2c libraries
import bioframe

import cooler
import cooltools

import glob
import os
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description='Plot log-binned contact vs distance for multiple samples.')
parser.add_argument('-i',required=True, metavar='mcool_dir',
                    help='directory containing .mcool files')
parser.add_argument('-o', required=True, metavar='outfile',
                    help='name of output figure file')
parser.add_argument('-a', required=False, action='store_true',
                    help='include all .mcool')

args = parser.parse_args()

# get filenames for all unfiltered mcools
if args.a:
    # use all mcool files in directory, regardless of name
    coolfiles = glob.glob(args.i + "/*.mcool")
else:
    # use for folders taht are direct output of distiller
    coolfiles = glob.glob(args.i + "/*.no_filter.1000.mcool")





# Use bioframe to fetch the genomic features from the UCSC.
mm9_chromsizes = bioframe.fetch_chromsizes('mm9')
mm9_cens = bioframe.fetch_centromeres('mm9')
# create a view with chromosome arms using chromosome sizes and definition of centromeres
mm9_arms = bioframe.make_chromarms(mm9_chromsizes, mm9_cens)

# get rid of chromarm in name. bc mouse chromosomes acrocentric
mm9_arms['name'] = mm9_arms['name'].apply(lambda x: x.rstrip('_p'))
mm9_arms = mm9_arms[mm9_arms['name']!='chrM']

# function for calculating ocntact v distance for one cool
def calculate_contacts_vs_distance(coolfile, resolution=1000):
    clr = cooler.Cooler(coolfile)

    # cvd == contacts-vs-distance
    cvd_smooth_agg = cooltools.expected_cis(
        clr=clr,
        view_df=mm9_arms,
        smooth=True,
        aggregate_smoothed=True,
        nproc=1 #if you do not have multiple cores available, set to 1
    )

    # Just take a single value for each genomic separation
    cvd_smooth_agg['s_bp'] = cvd_smooth_agg['dist']* resolution
    cvd_smooth_agg['balanced.avg.smoothed.agg'].loc[cvd_smooth_agg['dist'] < 2] = np.nan
    cvd_merged = cvd_smooth_agg.drop_duplicates(subset=['dist'])[['s_bp', 'balanced.avg.smoothed.agg']]
    
    return cvd_merged



df = pd.DataFrame()

# loop through each cool file
resolution=1000

for mcool in coolfiles:
    coolfile = mcool + '::/resolutions/'+str(resolution)
    print('Analyzing %s' % coolfile)
    temp = calculate_contacts_vs_distance(coolfile)
    temp['sample'] = os.path.basename(mcool).split('.')[0]
    df = df.append(temp)


# plot each sample with a different line    
cvd_plot = sns.lineplot(data=df.sort_values(by='sample'), 
                        x='s_bp', y='balanced.avg.smoothed.agg', alpha=0.8, hue='sample')
cvd_plot.set(xscale="log", yscale="log", 
    xlabel='genomic distance (bp)', ylabel='contact probability')
cvd_plot.minorticks_off()
cvd_plot.set_aspect(0.8)
sns.despine()

#save figure
plt.savefig(args.o, bbox_inches='tight')