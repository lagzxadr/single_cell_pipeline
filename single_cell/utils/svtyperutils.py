#Douglas Abrams
#10/14/19

import sys as sys
import numpy as np
import math as m
import pandas as pd

'''
takes in a dict and converts
to string
@info = dict
#from https://stackoverflow.com/questions/49693464/converting-dictionary-into-string/49694114
'''
def info_tostr(info):
    if type(info) == dict:
        return ";".join(("{}={}".format(*i) for i in info.items()))
    else:
        return info

'''
create an alt section for an output vcf
@param strand = string representing the strand senses
@param pos = pos of alt
@param chrom = chrom of alt
@param ref = reference sequence
'''
def makeAlt(strand, pos, chrom, ref):
    insert = str(chrom)+":"+str(pos)
    ref = str(ref)
    if strand == "--":
        return "["+insert+"["+ref
    elif strand == "++":
        return ref+"]"+insert+"]"
    elif strand == "-+":
        return "]"+insert+"]"+ref
    elif strand == "+-":
        return ref+"["+insert+"["
    else: return 0

'''
inserts row into pandas obj at
a specified index
@param data: pandas obj
@param new_row: row of data
    to be added
@param index: index at which to
    add row
'''
def addRow(data, new_row, index):

    before = data.iloc[:index]
    after = data.iloc[index:]

    with_new = before.append(new_row)

    return pd.concat([before, new_row,  after]).reset_index(drop = True)


'''
creates the info section for
the outputvcf from input destruct
csv
@param svtype: type of structural
    variation; see dict "to_svtype"
    in csv2vcf() for options
@param csv: input csv pandas obj
@param index: current index
'''
def expandINFOsection(svtype, loaderobj, index, variant_caller):

    start = loaderobj["POS"][index]
    end = loaderobj["POS2"][index]
    
    start_ci = "0,200"
    end_ci = "0,200"

    info = {'SVTYPE' : svtype,
            'CIPOS' : start_ci,
            'CIPOS95' : start_ci,
            'CIEND' : end_ci,
            'CIEND95' : end_ci }

    if svtype is "BND":
        info["MATEID"] = str(index)+"_2"

    else:
      #  info["SVLEN"] = end - start
        info["END"] = end
    return info


'''
parses a VCF file with svtyper annotation
to a csv file, or a csv file with svtyper
annotation to a vcf file
@param vcf: input vcf file/output filename
@param csv: input csv file/ output filename
'''
def extractSvtyperInfo(df):
    nRecords = len(df.index)
    data = []
    headers = None
    for i in range(nRecords):
        headers = df.iloc[i][-2].split(":")
        data.append(df.iloc[i][-1].split(":"))
    out = pd.DataFrame(data, columns = headers)
    return out

'''
for destruct csv*
adds partner entries for translocation calls
in output vcf
@param new_data: pandas obj representing output vcf
@param csv: pandas obj representing input destruct csv
@param nRecords: length on pandas objs
'''
def addBND_mates(new_data, loaderobj, nRecords, variant_caller):
    #loop through everything
    row_num = 0
    old_data_index = 0

    while row_num < nRecords:

        if new_data["INFO"][row_num]["SVTYPE"] == "BND":

            chrom = loaderobj["CHROM"][old_data_index]
            pos = loaderobj["POS"][old_data_index]
            id = str(row_num) + "_2"
            ref = new_data["REF"][row_num]
            strand =loaderobj["STRAND"][old_data_index]
            alt = makeAlt(strand, pos, chrom, ref)
            qual = "."
            filter = new_data["FILTER"][row_num]
            end = new_data["POS"][row_num]
            mate_id = new_data["ID"][row_num]
            info = expandINFOsection("BND", loaderobj, old_data_index, variant_caller)
            info["MATEID"] = str(row_num) + "_1"

       
            BND_mate = pd.DataFrame({'#CHROM':chrom,'POS':pos,'ID':id,
                                     'REF':ref,'ALT':alt,'QUAL':qual,
                                     'FILTER':filter,
                                     'INFO': info_tostr(info)},
                                     index = [row_num],
                                     columns = ['#CHROM','POS','ID',
                                                'REF', 'ALT', 'QUAL',
                                                'FILTER',
                                                'INFO'])

            new_data = addRow(new_data, BND_mate, row_num+1)
  
            nRecords+=1
            row_num+=1
        row_num+=1
        old_data_index+=1
        
    return new_data
'''
takes an inputcsv with an svtyper
annotation and writes it out to outfile
in matrix form.
@param str annotation = svtyper annotation to write
@param pd.df | str inputcsv = inputcsv containing svtyper annotation,
    cell_ids ["cell_id"] and chromosomes ["CHROM"], can be passed as pandas df or
    filename
@param str outfile = filename for write out file
@param str fillter = value for empty matrix indxs
'''
def writeSvtyperAnnotation(annotation, inputcsv, outfile, filler = np.nan):
    if str(inputcsv):
        inputcsv = pd.read_csv(inputcsv, delimiter = ",")
    nRecords = len(inputcsv.index)/2
    cell_ids = np.unique(inputcsv["cell_id"])
    
    positions = inputcsv["POS"]
    positions = positions[:nRecords]
    
    chroms = inputcsv["CHROM"]
    chroms = chroms[:nRecords]
    
    output = pd.DataFrame({"chromosome" : chroms, 
        "position" : positions}, columns = ["chromosome", "position"])

    for cell_id in cell_ids:
        svtyper_annotation = inputcsv.loc[inputcsv['cell_id'] == cell_id, annotation].values
        output[cell_id] = pd.Series(svtyper_annotation)
    output.to_csv(outfile,  index=False)




    
