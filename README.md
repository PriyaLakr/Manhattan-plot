# Manhattan-plot

Author: Priya Lakra


🐍 Python script for generating manhattan plot 🤓
Written in python 3

System requirements: No high performance computing required. Works well in a laptop with even 4GB RAM. 

Input -> File containing SNPs information such as their chromosomal location, p-values, etc. Currently acceptable formats include .csv, .tsv, .txt 

Output -> A manhattan plot in .png format with colored and black-white format 

What is a manhattan plot? 🤔 https://en.wikipedia.org/wiki/Manhattan_plot 

USAGE: 

    $ python pl_manhattan.py [options]
    usage: pl_manhattan.py [—-help] [--path path_to_input files] [--in_file input_file] 
                                      [--chromosome how_chromosome_is_denoted_in_input_file] 
                                      [--pval p_value] 
                                      [--position how_chromosomal_position_is_denoted_in_input_file] 
                                      [--col colored graph]  
    

Example:

Data: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013658/
Reference: Donkel SJ, Portilla Fernandez E, Ahmad S, Rivadeneira F, Van Rooij FJ, Ikram MA, Leebeek FW,  De Maat M, Ghanbari M (2021) Common and rare variants genetic association analysis of circulating neutrophil extracellular traps (NETs). Frontiers in Immunology. doi: 10.3389/fimmu.2021.615527.

    python pl_manhattan.py --path /Users/priyalakra/Desktop/Practice --in_file GCST90013658_buildGRCh37.tsv  --chromosome chromosome --pval p_value --position base_pair_location --col y

Output: GCST90013658_buildGRCh37_plot.png 


