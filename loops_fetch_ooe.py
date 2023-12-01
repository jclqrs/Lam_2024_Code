# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
from cooltools.api import snipping
import multiprocessing

# import libraries for biological data analysis
import cooler
import bioframe

import cooltools # cooltools v0.5.0
import cooltools.expected
from cooltools.api import snipping

# Set up parallelization
import multiprocess

import pybedtools
from pybedtools import BedTool
import argparse
import sys

parser = argparse.ArgumentParser(description='Fetch obs/exp_donut values for multiple contact maps.')
parser.add_argument('-c',required=True, metavar='config',
                    help='config file with name, map, and expected columns')
parser.add_argument('-l', required=True, metavar='looplist',
                    help='csv with unique combined_loop_id for each loop')
parser.add_argument('-r', required=True, metavar='resolution',
                    help='map resolution', type=int)
parser.add_argument('-p', required=True, metavar='nprocessors',
                    help='number of cores', default=1, type=int)

args = parser.parse_args()


ncpu = args.p
combined_loop_list = pd.read_csv(args.l)
mapres = args.r
config = args.c 
out = os.path.basename(args.c).split('.')[0]


# make sure each loop has a unique identifier
combined_loop_list['combined_loop_id'] = combined_loop_list.index

if 'combined_loop_id' not in combined_loop_list.columns.to_list():
    combined_loop_list['combined_loop_id'] = combined_loop_list.index
elif len(combined_loop_list.combined_loop_id.unique()) != len(combined_loop_list):
    sys.exit('each loop must have a unique combined_loop_id') 



# manually create list because bioframe is glitchy...
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
ends = [197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172, 129993255, 121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430, 166650296]

mm9_regions = pd.DataFrame()
mm9_regions['chrom'] = chroms
mm9_regions['start'] = 0
mm9_regions['end'] = ends 
mm9_regions['name'] = chroms


bedpecols = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']

def compare_looplists(loopsC, preC, loopsT, preT):
    """
    look at overlap of loops between two loop lists
    loopC, preC = output from dots function, name you want to give this loop list
    loopT, preT = output from dots function, name you want to give second loop list
    
    
    output: dataframe with merged list of loops and annotation of which maps they were called in
    """
    
    # read loop lists
    cols = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    cols2 = ['chrom1b', 'start1b', 'end1b', 'chrom2b', 'start2b', 'end2b']
    
    if type(loopsC) == str:
        dfC = pd.read_table(loopsC)[cols]
    else:
        dfC = loopsC[cols]
        
    if type(loopsT) == str: 
        dfT = pd.read_table(loopsT)[cols]
    else:
        dfT = loopsT[cols]

    # make bedtool for each loop list
    loopsC_bt = BedTool.from_dataframe(dfC)
    loopsT_bt = BedTool.from_dataframe(dfT)

    # intersect!
    # unique T
    uniqueT = loopsT_bt.pair_to_pair(loopsC_bt, 
                                       **{'type': 'notboth',
                                          'slop': 1,
                                          'is': True}).to_dataframe(names=cols)
    uniqueT['loopcall_map'] = preT

    #if shared, take loop peak at the control location
    shared = loopsC_bt.pair_to_pair(loopsT_bt, 
                                       **{'type': 'both',
                                          'slop': 1,
                                          'is': True}).to_dataframe(names=cols + cols2+ ['a', 'b', 'c', 'd'])
    shared = shared[cols]
    shared['loopcall_map'] = 'both'
    shared = shared.drop_duplicates()

    # no T loop called... 
    uniqueC = loopsC_bt.pair_to_pair(loopsT_bt, 
                                       **{'type': 'notboth',
                                          'slop': 1,
                                          'is': True}).to_dataframe(names=cols)
    uniqueC['loopcall_map'] = preC

    # merge C and T list
    combined_loop_list = uniqueT.append(shared).append(uniqueC)

    # label with unique loop id
    combined_loop_list['combined_loop_id'] = range(len(combined_loop_list))

    keep_cols = ['chrom1', 'start1', 'end1',
                'chrom2', 'start2', 'end2', 
                 'combined_loop_id', 'loopcall_map']

    combined_loop_list = combined_loop_list[keep_cols]

    return combined_loop_list


change_colors = {'weakened': 'firebrick',
                'strengthened': 'steelblue', 
                'unchanged': 'gray'}

def annotate_change(log2fc, t=1.2):
    if log2fc < np.log2(1/t):
        ret = 'weakened'
    elif log2fc > np.log2(t):
        ret = 'strengthened'
    else:
        ret = 'unchanged'
    return ret

def generate_ucsc_string(combined_loop_list):
    cols = combined_loop_list.columns
    udf = pd.DataFrame()
    udf['c'] = combined_loop_list['chrom1'].apply(lambda x: x+":")
    udf['s'] = combined_loop_list['start1'].apply(lambda x: "{:,}".format(x-500000)+'-')
    udf['e'] = combined_loop_list['start2'].apply(lambda x: "{:,}".format(x+500000))
    combined_loop_list['ucsc_window'] = udf.apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    return combined_loop_list.ucsc_window

def plot_annotated_matrix(mtx):
    plt.imshow(
        mtx,
        vmax = 4,
        vmin = 1,
        cmap='coolwarm')

    plt.colorbar(label = 'mean obs/exp')
    ticks_pixels = np.linspace(0, flank*2//resolution,5)
    ticks_kbp = ((ticks_pixels-ticks_pixels[-1]/2)*resolution//1000).astype(int)
    plt.xticks(ticks_pixels, ticks_kbp)
    plt.yticks(ticks_pixels, ticks_kbp)
    plt.xlabel('relative position, kbp')
    plt.ylabel('relative position, kbp')

    # Loop over data dimensions and create text annotations.
    for i in range(np.shape(mtx)[0]):
        for j in range(np.shape(mtx)[1]):
            text = plt.text(j, i, '%.2f' % (mtx[i, j]),
                           ha="center", va="center", color="w")


    plt.show()



def expand_coord(coord, res=10000, sort=0):
    # make array of coordinates surrounding + including peak
    shifts = np.array([-1,0,1, -1,0,1, -1,0,1]) *res
    expanded_pixels = np.repeat(coord, 9) + shifts
    if sort:
        expanded_pixels = np.sort(expanded_pixels)
        
    return expanded_pixels.astype(int)
 

def expand_peak_pixel(loop_series, res=10000):
    shifts = np.array([-1,0,1]) *res
    anc1 = loop_series[['start1', 'end1']] 
    anc2 = loop_series[['start2', 'end2']]

    expand1 = pd.DataFrame()
    expand2 = pd.DataFrame()

    for shift in list(shifts)*3:
        expand1 = expand1.append(anc1 + shift)

    for shift in list(shifts)*3:
        expand2 = expand2.append(anc2 + shift)

    expand1['start2'] = expand2.start2.sort_values().values
    expand1['end2'] = expand2.end2.sort_values().values

    # add chrom
    expand1 = expand1.astype(int)
    expand1['chrom1'] = loop_series['chrom1']
    expand1['chrom2'] = loop_series['chrom2']

    # reorder columns
    expand1 = expand1[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]
    
    return expand1


## To Do: also count nan values
# def count_nans(loops, clr, expected, mm9_regions, resolution=10000, flank=50000):
#     stack = cooltools.pileup(clr, 
#                     loops, 
#                     expected_df=expected, 
#                     flank=flank,
#                     view_df=mm9_regions, nproc=ncpu)
    
#     # mark nan values
#     nan_values = np.isnan(stack)
    
#     # count for each loop + flank region
#     nan_count = np.sum(nan_values, axis=(0,1))
    
#     return nan_count

## To Do: also count nan values
def calc_obs(loops, clr,mm9_regions, resolution=10000, flank=10000):
    stack = cooltools.pileup(clr, 
                    loops, 
                    flank=flank,
                    view_df=mm9_regions, nproc=ncpu)
    
    # count for each loop + flank region
    obs_count = np.sum(stack, axis=(0,1))
    
    return obs_count

def calc_looe_donuts(loops, clr, exp, mm9_regions, resolution):

    flank = resolution * 5
    
    # generate round kernel
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
    
    loops_expand = make_expanded_df(loops, resolution=resolution)
    loops_expand['region'] = loops_expand.chrom1
    

    # Make pileup for each heatmap
    # cooltools pileup function will determine center point of each anchor
    # and then add/subract flank value to the center
    stack = cooltools.pileup(clr, 
                    loops_expand[[ 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']], 
                    expected_df=exp, 
                    flank=flank,view_df=mm9_regions, nproc=ncpu)
    
    # grab peak pixel
    loops_expand['ooe'] = stack[5,5,:]
    
    
    
    # grab kernel_o/e
    
    snipper = snipping.CoolerSnipper(
    clr,
    view_df=mm9_regions,
    cooler_opts={"balance": "weight"},
    min_diag=2
    )


    snipper2 = snipping.ExpectedSnipper(
        clr, exp,
        view_df=mm9_regions,
        min_diag=2
    )


    features_df = loops_expand.copy()
    features_df = snipping.expand_align_features(
        features_df, flank, clr.binsize, format='bedpe'
    )


    features_df["lo1"] = (features_df["start1"] / clr.binsize).astype(int)
    features_df["hi1"] = (features_df["end1"] / clr.binsize).astype(int)
    features_df["lo2"] = (features_df["start2"] / clr.binsize).astype(int)
    features_df["hi2"] = (features_df["end2"] / clr.binsize).astype(int)


    # Find region offsets and then subtract them from the feature extents

    region_offsets = mm9_regions[["chrom", "start", "end"]].apply(clr.offset, axis=1)
    region_offsets_dict = dict(zip(mm9_regions["name"].values, region_offsets))

    features_df["region_offset"] = features_df["region"].replace(region_offsets_dict)



    features_df[["lo1", "hi1"]] = (
        features_df[["lo1", "hi1"]]
        .subtract(features_df["region_offset"].fillna(0), axis=0,)
        .astype(int)
    )
    features_df[["lo2", "hi2"]] = (
        features_df[["lo2", "hi2"]]
        .subtract(features_df["region_offset"].fillna(0), axis=0,)
        .astype(int)
    )
    

    pool = multiprocessing.Pool(ncpu)
    mymap = pool.map

    stackO = snipping.pileup_legacy(features_df, snipper.select, snipper.snip, map=mymap)
    stackE = snipping.pileup_legacy(features_df, snipper2.select, snipper2.snip, map=mymap)

      
    stacked_donut = np.dstack([mask] * np.shape(stackO)[2])
    stackO_donut = stackO * stacked_donut
    stackE_donut = stackE * stacked_donut

    loops_expand['obs'] = stackO[5,5,:]
    
    loops_expand['rounddonut_oe'] = np.nansum(stackO_donut, axis=(0,1)) / np.nansum(stackE_donut, axis=(0,1))

    loops_expand['ooe_donut'] = loops_expand['ooe'] / loops_expand['rounddonut_oe']
    
    return loops_expand



def make_expanded_df(loop_df, resolution):
    test = loop_df.copy()
    test['a'] = test['start1'].apply(expand_coord, res=resolution)
    test['b'] = test['end1'].apply(expand_coord, res=resolution)
    test['c'] = test['start2'].apply(expand_coord, sort=1, res=resolution)
    test['d'] = test['end2'].apply(expand_coord, sort=1, res=resolution)
    
    test = test.explode(list('abcd'))
    
    test = test[['chrom1', 'a', 'b', 'chrom2', 'c', 'd', 'combined_loop_id']]
        
    test.columns = [ 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'combined_loop_id']
    
    test[['start1', 'end1', 'start2', 'end2']] = test[['start1', 'end1', 'start2', 'end2']].astype(int)
    
    return test




def calc_looe_donuts_PEAK(loops, clr, exp, mm9_regions, resolution):

    flank = resolution *5
    
    # generate round kernel
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
    
    
    # since we're only looking at the peak, don't expand pixels
    loops_expand = loops.copy()
    loops_expand['region'] = loops_expand.chrom1
    

    # Make pileup for each heatmap
    # cooltools pileup function will determine center point of each anchor
    # and then add/subract flank value to the center


    stack = cooltools.pileup(clr, 
                    loops_expand[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']], 
                    expected_df=exp, 
                    flank=flank,view_df=mm9_regions, nproc=ncpu)
    
        
    # # mark nan values
    # nan_values = np.isnan(stack)
    
    # # count for each loop + flank region
    # nan_count = np.sum(nan_values, axis=(0,1))
    

    # switch to counting zero values... seems like a better measure of sparsity
    nan_values = (np.nan_to_num(stack)==0)
    
    # count for each loop + flank region
    nan_count = np.sum(nan_values, axis=(0,1))

    # grab peak pixel
    loops_expand['ooe'] = stack[5,5,:]
    
    # grab kernel_o/e
    
    snipper = snipping.CoolerSnipper(
    clr,
    view_df=mm9_regions,
    cooler_opts={"balance": "weight"},
    min_diag=2,
    )


    snipper2 = snipping.ExpectedSnipper(
        clr, exp,
        view_df=mm9_regions,
        min_diag=2,
    )


    features_df = loops_expand.copy()
    features_df = snipping.expand_align_features(
        features_df, flank, clr.binsize, format='bedpe'
    )


    features_df["lo1"] = (features_df["start1"] / clr.binsize).astype(int)
    features_df["hi1"] = (features_df["end1"] / clr.binsize).astype(int)
    features_df["lo2"] = (features_df["start2"] / clr.binsize).astype(int)
    features_df["hi2"] = (features_df["end2"] / clr.binsize).astype(int)



    # Find region offsets and then subtract them from the feature extents

    region_offsets = mm9_regions[["chrom", "start", "end"]].apply(clr.offset, axis=1)
    region_offsets_dict = dict(zip(mm9_regions["name"].values, region_offsets))

    features_df["region_offset"] = features_df["region"].replace(region_offsets_dict)



    features_df[["lo1", "hi1"]] = (
        features_df[["lo1", "hi1"]]
        .subtract(features_df["region_offset"].fillna(0), axis=0,)
        .astype(int)
    )
    features_df[["lo2", "hi2"]] = (
        features_df[["lo2", "hi2"]]
        .subtract(features_df["region_offset"].fillna(0), axis=0,)
        .astype(int)
    )
    

    pool = multiprocessing.Pool(ncpu)
    mymap = pool.map

    stackO = snipping.pileup_legacy(features_df, snipper.select, snipper.snip, map=mymap)
    stackE = snipping.pileup_legacy(features_df, snipper2.select, snipper2.snip, map=mymap)

      
    stacked_donut = np.dstack([mask] * np.shape(stackO)[2])
    stackO_donut = stackO * stacked_donut
    stackE_donut = stackE * stacked_donut
    
    loops_expand['rounddonut_oe'] = np.nansum(stackO_donut, axis=(0,1)) / np.nansum(stackE_donut, axis=(0,1))

    loops_expand['ooe_donut'] = loops_expand['ooe'] / loops_expand['rounddonut_oe']
    
    loops_expand['nancount'] = nan_count

    loops_expand['obs'] = stackO[5,5,:]

    
    return loops_expand




df = pd.read_csv(config)

for contactmap, coolpath, exp in zip(df['name'], df['map'], df['expected']):
    print('Analyzing ' + contactmap)
    # Read in clrs and expected table
    clrC = cooler.Cooler(coolpath + "::/resolutions/" + str(mapres))
    expectedC = pd.read_table(exp)
    

    # calc  ooe_donut based on peak pixel

    loops_expandC = calc_looe_donuts_PEAK(combined_loop_list, clrC, expectedC, 
                                                                mm9_regions, resolution=mapres)
    cols = ['combined_loop_id', 'obs', 'ooe', 'ooe_donut', 'nancount']
    loops_expandC = loops_expandC[cols]
    loops_expandC.columns = ['combined_loop_id', 'peak_obs_' + contactmap, 'peak_ooe_' + contactmap, 'peak_ooe_donut_' + contactmap, 'nancount_' + contactmap]

    combined_loop_list = combined_loop_list.merge(loops_expandC, on='combined_loop_id', how='left')

    # calc based on expanded loop (9 pixel avg)
    # calc  ooe_donut for untreated sample

    loops_expandC = calc_looe_donuts(combined_loop_list, clrC, expectedC, 
                                                                mm9_regions, resolution=mapres)
    loops_expandC = loops_expandC.groupby('combined_loop_id')['ooe_donut'].mean().reset_index()
    loops_expandC.columns = ['combined_loop_id', 'mean_ooe_donut_' + contactmap]
    combined_loop_list = combined_loop_list.merge(loops_expandC, on='combined_loop_id', how='left')


combined_loop_list.to_csv(out + '_ooe_table.csv', index=None)