#run: pytest test_manhattan.py
#importing 
import pytest
import pandas as pd
import numpy as np
from pl_manhattan import data_processing
from pl_manhattan import prepare_file
from pl_manhattan import edit


def data_processing(file,chromosome):
    
    #create a new column for chromosomes data where they are stored as integers 
    #replace "X" chromosome as an integer number 
    file["new_chr"] = file[chromosome].replace(['X', 'Y'], ['23', '24'])
    file["new_chr"] = file["new_chr"].astype(int)
    
    return file
    
    
file_X = pd.DataFrame({ "values": ['A','T','G'], "chrom": [1,2,'X'] })
res_X = pd.DataFrame({ "values": ['A','T','G'], "chrom": [1,2,'X'], "new_chr": [1,2,23] })

file_Y = pd.DataFrame({ "values": ['A','T','G'], "chrom": [1,2,'Y'] })
res_Y = pd.DataFrame({ "values": ['A','T','G'], "chrom": [1,2,'Y'], "new_chr": [1,2,24] })

file_XY = pd.DataFrame({ "values": ['A','T','G','C'], "chrom": [2,1,'Y','X'], "p": [0,0.2,1000,0.00004] })
res_XY = pd.DataFrame({ "values": ['A','T','G','C'], "chrom": [2,1,'Y','X'], "p": [0,0.2,1000,0.00004], "new_chr": [2,1,24,23] })
res_log = pd.DataFrame({ "values": ['T','A','C','G'], "chrom": [1,2,'X','Y'], "p": [0.2,0,0.00004,1000], "new_chr": [1,2,23,24], "index": [2,1,4,3], "-log10(p)": [0.6989700043360187, 0, 4.3979400086720375, -3.0] })

def test_data_processing_X():
    result_X = data_processing(file=file_X, chromosome='chrom') 
    pd.testing.assert_frame_equal(res_X, result_X)
    
    
def test_data_processing_Y():
    result_Y = data_processing(file=file_Y, chromosome='chrom') 
    pd.testing.assert_frame_equal(res_Y, result_Y)

# the below test function throws error; 
# FAILED test_manhattan.py::test_data_processing - 
# ValueError: The truth value of a DataFrame is ambiguous. Use a.empty, a.bool(), a.item(), a....
# def test_data_processing(self):
    # assert data_processing(file=file2, chromosome='chrom') == result2
 
    
def test_data_processing_XY():
    result_XY = data_processing(file=file_XY, chromosome='chrom') 
    pd.testing.assert_frame_equal(res_XY, result_XY)

   
def test_prepare_file():
    result_log = prepare_file(file=res_XY, pval='p', chromosome='new_chr') 
    # dropping index as they don't need to be equal in this case
    pd.testing.assert_frame_equal(res_log.reset_index(drop=True), result_log.reset_index(drop=True)) 

    
def edit(file_new,position,chromosome):

    # add new position for plotting across x-axis
    positions=file_new[[chromosome, position]]  
    new_pos = []
    add=0
    for chro,posi in positions.groupby(chromosome):
        new_pos.append(posi[[position]]+add)
        add+=posi[position].max() # maximum position within each chromosome
        
    # append new positions to the data frame     
    file_new['new_pos'] = pd.concat(new_pos)
   
    return file_new
    
def test_edit():
    result_edit = edit(file=res_log, chromosome='chrom') 
  #  pd.testing.assert_frame_equal(res_Y, result_Y)
    
    
    