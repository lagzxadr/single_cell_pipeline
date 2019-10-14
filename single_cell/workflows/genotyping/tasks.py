#Douglas Abrams
#10/14/19

import pandas as pd
import single_cell.utils.helpers as helpers 
import os as os
import pypeliner
import vcf as vf
import csv as cv
import single_cell.utils.vcfutils as vcfutils
import single_cell.utils.svtyperutils as svtyperutils
import vcaller_loader as loader

'''
calls svtyper-sso on input
bam and vcf to perform genotyping.
'''
def genotype(input_bam, reference, input_vcf, 
        output_vcf, output_csv, tempdir, 
        docker_image = None):

        helpers.makedirs(tempdir)

        cmd = ['svtyper-sso', 
                '--input_vcf', input_vcf,
                '--bam', input_bam,
                '--ref_fasta', reference,
                '-o', output_vcf]
        base_data = vcfutils.parseVCF(input_vcf, None, return_pandas = True)
        
        pypeliner.commandline.execute(*cmd, docker_image = docker_image)

        base_data = vcfutils.parseVCF(output_vcf, None, return_pandas = True)

        svtype_annotations = svtyperutils.extractSvtyperInfo(base_data)

        base_data = base_data.iloc[:, :-2] #assumes svtyper info in last 2 cols

        output = pd.concat([base_data, svtype_annotations], axis = 1)

        output.to_csv(output_csv, index = False, encoding = "utf-8-sig")

'''
merges input csv files
into one csv
'''
def merge_csvs(input_csvs, merged_csv):
        helpers.makedirs(tempdir)
        vals = []

        if type(input_csvs) == dict:
                for key in input_csvs.keys():
                        vals.append(input_csvs[key])

        else:
                vals = input_csvs
                pd_concat = pd.concat([pd.read_csv(f) for f in vals],
                                                         sort = False)


        csv_output = open(merged_csv, mode ="w+")

        csv_writer = cv.writer(csv_output, delimiter = ',',
        quotechar = '"',
        quoting = cv.QUOTE_MINIMAL)

        pd_concat.to_csv(merged_csv,index=False,  encoding = "utf-8-sig")


'''
convert input variant call file csv
[either lumpy or destruct]
to vcf that svtyper can take
@param input -> input variant call CSV
@param vcf -> output path for vcf
@param tempdir -> temp dir in which to output temps
@param caller -> variant caller, either lumpy or destruct
'''
def varcalls_to_svtyper_input(input, vcf, tempdir, caller):

        helpers.makedirs(tempdir)       
  
        csv = loader.VarCall_Loader(input, caller)
        nRecords = len(csv.data.index)

        CHROM = csv["CHROM"]

        POS = csv["POS"]

        REF = ["N"]*nRecords

        QUAL = ["."] * nRecords

        FILTER = ["Pass"] * nRecords

        SVTYPE = csv["TYPE"]

        ID   = [None] * nRecords
        ALT  = [None] * nRecords
        INFO = [None] * nRecords
        
        for i in range(nRecords):      
                if SVTYPE[i] == "BND":
                        ID[i] = str(i)+"_1"
                        ALT[i] = svtyperutils.makeAlt(csv["STRAND"][i], 
                                                        POS[i], CHROM[i], 
                                                        REF[i])
                else:
                        ID[i] = str(i) 
                        ALT[i] = "<" + str(SVTYPE[i]) + ">"

                INFO[i] = svtyperutils.expandINFOsection(SVTYPE[i],  
                                                        csv, i, caller)
        
                new_data = pd.DataFrame({'#CHROM':CHROM, 'POS': POS, "ID":ID, 
                                         "REF": REF, 'ALT': ALT, 'QUAL': QUAL,
                                         'FILTER':FILTER, 'INFO':INFO},
                                columns = ['#CHROM', 'POS', 'ID', 'REF',
                                           'ALT', 'QUAL', 'FILTER', 'INFO'])

        new_data = svtyperutils.addBND_mates(new_data, csv, nRecords, caller)
        
        new_data["INFO"] = [svtyperutils.info_tostr(info) for info in new_data["INFO"]]


        new_data.to_csv(vcf, sep = "\t", index=False)



'''
writes the annotations contained in the below
annotations list to files, each to their own
@param csv -> csv file containg annotations
as features
@param output_dirs -> output directories for 
annotation files
'''
def writeSvtyperAnnotations(csv, output_paths,  tempdir):

        helpers.makedirs(tempdir)    

        annotations = ["AO", "AP", "AS",
                        "ASC","DP", "GQ",
                        "QA", "QR", "RO",
                        "RP", "RS", "SQ",
                        "GL", "AB"]
        
        

        for i in range(len(output_paths)):

                output_path = output_paths[i]
                annotation = annotations[i]

                svtyperutils.writeSvtyperAnnotation(annotation, 
                                                        csv,
                                                        output_path)



