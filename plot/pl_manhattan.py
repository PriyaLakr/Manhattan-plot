# importing 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import os, sys
import argparse 
from time import strftime 


def open_infile(in_file):
    
    if in_file.endswith(".csv"):
    	file_name=pd.read_csv(in_file,sep=',')
    elif in_file.endswith(".tsv") or in_file.endswith(".txt"):
        file_name=pd.read_csv(in_file,sep='\t')
    else:
    	# print the concern and exit the program 
        print("Upload a properly formatted file!")
        sys.exit(0) # or sys.exit(1)
            
    return file_name


def data_processing(file,chromosome):
    
    #create a new column for chromosomes data where they are stored as integers 
    #replace "X" chromosome as an integer number 
    file["new_chr"] = file[chromosome].replace(['X', 'Y'], ['23', '24'])
    file["new_chr"] = file["new_chr"].astype(int)
    
    return file


def prepare_file(file,pval,chromosome):
    
    # add an index column 
    file["index"] = [i for i in range(1,len(file[[chromosome]])+1)]
    
    # perform p value log transformation
    file["-log10(p)"] = np.where(file[[pval]] > 0, -np.log10(file[[pval]]), 0)
    
    file_new = file.sort_values(chromosome)
    
    return file_new
    
    
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


def manhattanplot(snpfile,col,chromosome):
      
    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_subplot(111)
    ax_label = []
    ax_pos = []
    if col == 'y':
        colors = itertools.cycle(['red','blue','green','gray','black','orange','beige','gold','coral','magenta'])
    else:
        colors = itertools.cycle(['gray','black'])

    for (chrname,detail) in list(snpfile.groupby(chromosome)): 
        detail.plot(kind = 'scatter', x = "new_pos", y = "-log10(p)", color = next(colors), ax=ax, s = 4);
        ax_label.append(chrname)
        ax_pos.append((detail['new_pos'].iloc[-1] + detail['new_pos'].iloc[0])/2)
    
    ax.set_xticks(ax_pos)
    ax.set_xticklabels(ax_label,rotation='vertical')
    ax.set_xlabel("Chromosomes")
    ax.set_ylabel("$-log_{10}$(p) values")
    ax.set_ylim(0,snpfile["-log10(p)"].max()+5)
    plt.axhline(y = -np.log10(1e-05), color = "magenta", linewidth = 0.5)
    plt.axhline(y = -np.log10(5e-08), color = "blue", linewidth = 0.5)

    plt.savefig(f'{filename}_plot.png')


if __name__ == "__main__":

	# initialize your parser
    parser = argparse.ArgumentParser(description = "Script for designing a manhattan plot")   

	# parse the arguments
    parser.add_argument("--path", type=str, help="Path of the directory where your SNP/GWAS files are stored")
    parser.add_argument("--infile", type=str, help="Name of the file containing SNPs/GWAS data to be analysed")
    parser.add_argument("--chromosome",  type=str, help="How chromosomes are indicated in your input file")
    parser.add_argument("--pval",  type=str, help="How p-values are indicated in your input file")
    parser.add_argument("--position", type=str, help="How chromosomal positions are indicated in your input file")
    parser.add_argument("--have_X", type=str, help="Does your data contain X chromosome? If yes, Write 'y' ")
    parser.add_argument("--col", type=str, help="Do you want a colored graph? If yes, Write 'y' ")

    args = parser.parse_args()
    have_X = args.have_X
    
    start_time = strftime("%m-%d-%Y %H:%M:%S")
    sys.stdout.write("Program started at " + start_time + '\n') # can also use print() here
    os.chdir(args.path)

    # call functions
    file = open_infile(args.infile) 

    if have_X == 'y':
        proc_file = data_processing(file,args.chromosome)
        filename = args.infile.rsplit(".",1)[0]
        new_path = f'{args.path}/{filename}_result'
        try:
        	os.makedirs(new_path)
        except OSError:
        	print('directory exists')
        os.chdir(new_path)
        
        proc_snpinfile = prepare_file(proc_file,args.pval,chromosome="new_chr") 
        proc_snpfile = edit(proc_snpinfile,args.position,chromosome="new_chr")
        manhattanplot(proc_snpfile,args.col,chromosome="new_chr") 

    else:
    # making a result folder to store results
        filename = args.infile.rsplit(".",1)[0]
        new_path = f'{args.path}/{filename}.result'
        try:
        	os.makedirs(new_path)
        except OSError:
        	print('dir exists')
        os.chdir(new_path)
		
        snpinfile = prepare_file(file,args.pval,args.chromosome) 
        snpfile = edit(proc_snpinfile,args.position,args.chromosome)
        manhattanplot(snpfile,args.col,args.chromosome)
        
    end_time = strftime("%m-%d-%Y %H:%M:%S")
    sys.stdout.write("Program ended at " + end_time + '\n')
    sys.exit(0)
