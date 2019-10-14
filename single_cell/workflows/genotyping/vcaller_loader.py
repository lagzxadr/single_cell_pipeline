import enum
import sys
import pandas as pd

svtyper_set = ["CHROM", "POS", "CHROM2", 
               "POS2"] #and "STRAND", "ID", "ALT"

lumpy_set = ["chrom1", "start1", "chrom2", 
             "start2"]

destruct_set = ["chromosome_1", "position_1", 
                "chromosome_2", 
                "position_2"]

lumpy_to_svtype = {"INTERCHROM":"BND",
             "DUPLICATION":"DUP",
             "DELETION":"DEL",
             "INSERTION":"INS",
             "INVERSION":"INV"}

destruct_to_svtype =  {"translocation":"BND",
                       "duplication":"DUP",
                       "deletion":"DEL",
                       "insertion":"INS",
                       "inversion":"INV"}
'''
abstracts out feature names for specific file types
uses "svtyper_set" usernames to pull type-specific username
when indexing. Overloads __getitem___
@param csv -> input csv, either lumpy or destruct
@param caller -> name of variant caller, either "lumpy" or "destruct
@param delim -> delim for csv read, auto set so that for 
    lumpy -> "," and for destruct -> "\t"
'''
class VarCall_Loader:

    def __init__(self, csv, caller, delim = None):
        if not delim:
            if caller is "lumpy": delim = ","
            if caller is "destruct": delim = "\t"

        self.data =  pd.read_csv(csv, delim)
        self.variant_caller = caller
        
        self.translator = dict(zip(svtyper_set, destruct_set))

        if caller == "lumpy":
            self.translator = dict(zip(svtyper_set, lumpy_set))

    def __getitem__(self, index):        
        if index == "STRAND" and self.variant_caller == "lumpy":
            return self.data["strands"]

        if index == "STRAND" and self.variant_caller == "destruct":
            return self.data["strand_1"].combine(self.data["strand_2"], 
                        lambda s1, s2: str(str(s1) + str(s2))) 

        if index == "TYPE" and self.variant_caller == "lumpy":  
            return  self.data["type"].apply(lambda type: lumpy_to_svtype[type])

        if index == "TYPE" and self.variant_caller == "destruct":   
            return  self.data["type"].apply(lambda type: destruct_to_svtype[type])

        else: return self.data[self.translator[index]]



# def test():
#     vcl = VarCall_Loader(sys.argv[1], "destruct")
#     print vcl["TYPE"]
    
# if __name__ == "__main__":
#     test()