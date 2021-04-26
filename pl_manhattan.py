import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import os
import argparse 

def open_infile(infile):
    
    if infile.endswith(".csv"):
    	file_name=pd.read_csv(infile,sep=',')
    elif infile.endswith(".tsv") or infile.endswith(".txt"):
        file_name=pd.read_csv(infile,sep='\t')
    else:
	print("Incorrect file format")            
    return file_name


def prepare_file(file,chromosome,pval,position):
    
    # add an index column 
    file["index"] = [i for i in range(1,len(file[[chromosome]])+1)]
    
    # perform p value log transformation
    file["-log10(p)"] = -np.log10(file[[pval]])

    file_new = file.sort_values(chromosome)
    
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


def manhattanplot(snpinfile):
    
    
    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax_label = []
    ax_pos = []
    colors = itertools.cycle(['red','blue','green','gray','black','orange','beige','gold','coral','magenta'])


    for (chrname,detail) in list(snpinfile.groupby(chromosome)): 
        detail.plot(kind = 'scatter', x = "new_pos", y = "-log10(p)", color = next(colors), ax=ax);
        ax_label.append(chrname)
        ax_pos.append((detail['new_pos'].iloc[-1] + detail['new_pos'].iloc[0])/2)
    
    ax.set_xticks(ax_pos)
    ax.set_xticklabels(ax_label,rotation='vertical')
    ax.set_xlabel("Chromosomes")
    ax.set_ylabel("$-log_{10}$(p) values")
    ax.set_ylim(0,snpinfile["-log10(p)"].max()+5)
    plt.axhline(y = -np.log10(1e-05), color = "magenta", linewidth = 0.5)
    plt.axhline(y = -np.log10(5e-08), color = "blue", linewidth = 0.5)

   # plt.show()

    plt.savefig(f'{filename}_plot.png')


if __name__ == "__main__":
	# initialize your parser
    parser = argparse.ArgumentParser(description = "This script design molecular inversion probes")   

	# parse the arguments
    parser.add_argument("--path", type=str, help="Path of the directory where SNP csv files are stored")
    parser.add_argument("--in_file", type=str, help="Name of the csv and/or tsv file containing SNPs/GWAS data to be analysed")
    parser.add_argument("--chromosome",  type=str, help="How chromosomes are indicated in csv file")
    parser.add_argument("--pval",  type=str, help="How p-values are indicated in csv file")
    parser.add_argument("--position", type=str, help="How chromosomal positions are indicated in csv file")

    args = parser.parse_args()
    path = args.path
    infile = args.in_file
    chromosome = args.chromosome
    pval = args.pval
    position = args.position

    os.chdir(path)

    # call functions
    file=open_infile(infile)

    # making a result folder to store results
    filename = infile.rsplit(".",1)[0]
    new_path = f'{path}/{filename}.result'
    os.makedirs(new_path)
    os.chdir(new_path)
		
    snpinfile=prepare_file(file,chromosome,pval,position) 
    manhattanplot(snpinfile) 






