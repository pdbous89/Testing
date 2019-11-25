#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 11:03:05 2019
@author: pbousounis
"""

import os
import numpy as np
import pandas as pd
import pybedtools
from pybedtools import BedTool

# get today's date as YYYY-MM-DD
today = datetime.today().strftime('%Y-%m-%d')

# Review the ClinVar bed file
#cv_bed_colnames = ['chr', 'start', 'end', 'name']
cv_bed_file = '/Users/pbousounis/Experiments/2019-10-29_hg19mod/2019-10-29-RefSeq-CV_recovery/2019-11-05_clinvar_path-l.bed'
cv_bed = pd.read_csv(cv_bed_file, sep='\t', header=None)
cv_bed.head()

# Review the RefSeq exons file
rs_exons_file = '/Users/pbousounis/Experiments/2019-10-29_hg19mod/2019-10-25_RefSeq_GFF3_download_validate/2019-11-05_GRCh37_latest_genomic.gff_bed/2019-11-05_GRCh37_latest_genomic.gff.bed'
rs_exons = pd.read_csv(rs_exons_file, sep='\t', header=None, low_memory=False)
rs_exons = rs_exons.drop([0])
rs_exons.head()

#rs_exons.to_csv(rs_exons_file, sep='\t', header=None, index=False)

# define the not-element-of function to get regions unique to the test_bed file
def bed_not(ref_bed_filepath, test_bed_filepath):
    
    """ Given two bed files, ref_bed and test_bed, perform bedtools intersect -c to return the number of 
    test_bed regions that DO overlap any regions in the ref_bed file. Returns the ref_bede file with a column of overlap counts"""
    
    cwd = os.getcwd()
    
    # specify the reference bed file
    ref_bedtool = BedTool(ref_bed_filepath)
    prfx_ref = ref_bed_filepath.split('/')[-1]
    prfx_ref = prfx_ref.split('.')[0]
    
    # specify the new ClinVar bed file
    test_bedtool = BedTool(test_bed_filepath)
    prfx_test = test_bed_filepath.split('/')[-1]
    prfx_test = prfx_test.split('.')[0]

    # specify name/path of output bed file
    out_bed_filepath = '{}/{}_IN_{}.bed'.format(cwd, prfx_test, prfx_ref)
    
    # run bedtools intersect to get all test_bed regions NOT found in ref_bed (-v option)
    ref_NotIn_test = test_bedtool.intersect(b=ref_bedtool, c=True)
    ref_NotIn_test.head()
    
    ref_NotIn_test.saveas(out_bed_filepath, trackline="track name='ClinVar loci IN TSO bed file'", )
    feat_count = ref_NotIn_test.count()
    
    # confirm file saved
    print('\nNumber of {} features NOT IN {} = {}\n'.format(prfx_test, prfx_ref, feat_count))
    
    if os.path.exists(out_bed_filepath):
        print('Success!\nFile saved to: \n{}.'.format(os.path.join(cwd, out_bed_filepath)))

    return(feat_count)
    
    
CvRs_bednot = bed_not(cv_bed_file, rs_exons_file)