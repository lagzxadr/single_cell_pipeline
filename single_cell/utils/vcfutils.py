'''
Created on Feb 27, 2018

@author: dgrewal
'''
import os
import logging
import numpy as np
import sys
import math as m
import random
import string
import os.path
import pandas as pd
from os import path
import helpers as helpers
import csv as cv
'''
parses a VCF file to a CSV file.
@param vcf: input vcf file
@param csv: output csv filename
@param return_pandas: if True, 
    return pandas obj, make no file
'''
def parseVCF(vcf, csv, return_pandas = False):
    data = []
    cols = []
    with open(vcf) as fp:
        for cnt, line in enumerate(fp):
            if line[0] != "#":
                data.append( line.strip().split("\t"))
            if line[0] == "#" and line[1] != "#":
                cols = line.strip().split("\t")
                if cols[0] == "#CHROM": cols[0] = "CHROM"
    vcf = pd.DataFrame(data,  columns = cols)
    if return_pandas:
        return vcf
    if not return_pandas:
        csv = pd.read_csv(csv, delimiter = "\t")
        new_data.to_csv(csv,  index=False, encoding = "utf-8-sig")

'''
csv to vcf
'''
def parseCSV(csv, vcf, return_pandas = False):
    data = pd.read_csv(csv, delimiter = ",")
    if outputPandas:
        return data
    else:
        data.to_csv(vcf, sep = "\t", index = False)
        
    return out

'''
adds key label supplementary
info in file dict to files
'''
def addSuppInfo(k, v, label):
    df = pd.read_csv(v, delimiter = ",")
    supp_col = [k]*len(df.index)
    df[label] = supp_col
    print df
    df.to_csv(v, sep = ",", index = False)

'''
merges input csv files
into one csv
'''

def merge_csvs(input_csvs, merged_csv, tempdir = None):
    if tempdir: helpers.makedirs(tempdir)

    vals = []

    csv_output = open(merged_csv, mode ="w+")

    csv_writer = cv.writer(csv_output, delimiter = ',',
                                           quotechar = '"',
                                           quoting = cv.QUOTE_MINIMAL)

    if type(input_csvs) is dict:

        for key in input_csvs.keys():
            addSuppInfo(key, input_csvs[key], "cell_id")
            vals.append(input_csvs[key])
            pd_concat = pd.concat([pd.read_csv(f) for f in vals], sort = False)            
    else:
        vals = input_csvs
        val1 = input_csvs[0]
        val2 = input_csvs[1]
        
        one = pd.read_csv(val1)
        chunks =  pd.read_csv(val2, chunksize = 50000)
        for chunk in pd.read_csv(val2, chunksize = 50000):
            one = pd.concat([one, chunk] , sort = False)    
        pd_concat = one



    pd_concat.to_csv(merged_csv,index=False,  encoding = "utf-8-sig")


#########################################
#DOUGLAS CODE ENDS
#########################################



def _get_header(infile):
    '''
    Extract header from the VCF file

    :param infile: input VCF file
    :return: header
    '''

    header = []
    for line in infile:
        if line.startswith('##'):
            header.append(line)
        elif line.startswith('#'):
            header.append(line)
            return header
        else:
            raise Exception('invalid header: missing #CHROM line')

    logging.getLogger("single_cell.helpers.vcfutils").warn(
        "One of the input files is empty"
    )
    return []

def concatenate_vcf(infiles, outfile):
    '''
    Concatenate VCF files

    :param infiles: dictionary of input VCF files to be concatenated
    :param outfile: output VCF file
    '''

    with open(outfile, 'w') as ofile:
        header = None

        for _,ifile in infiles.iteritems():

            if os.path.getsize(ifile) == 0:
                logging.getLogger("single_cell.helpers.vcfutils").warn(
                    'input file {} is empty'.format(ifile)
                )
                continue

            with open(ifile) as f:

                if not header:
                    header = _get_header(f)

                    for line in header:
                        ofile.write(line)
                else:
                    if not _get_header(f) == header:
                        logging.getLogger("single_cell.helpers.vcfutils").warn(
                            'merging vcf files with mismatching headers'
                        )

                for l in f:
                    print >> ofile, l,

