# Manhattan-plot

🐍 Python script for generating manhattan plot 🤓

Input -> csv file containing SNPs information such as their chromosomal location, p-values, etc. Currently acceptable formats include .csv, .tsv, .txt 
Output -> A manhattan plot in .png format with colored and black-white format 

What is a manhattan plot? 🤔 https://en.wikipedia.org/wiki/Manhattan_plot 

USAGE help: 
python pl_manhattan.py --help

Example:

Data: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013658/
Reference: Donkel SJ, Portilla Fernandez E, Ahmad S, Rivadeneira F, Van Rooij FJ, Ikram MA, Leebeek FW,  De Maat M, Ghanbari M (2021) Common and rare variants genetic association analysis of circulating neutrophil extracellular traps (NETs). Frontiers in Immunology. doi: 10.3389/fimmu.2021.615527.

python pl_manhattan.py --path /Users/priyalakra/Desktop/Practice --in_file GCST90013658_buildGRCh37.tsv  --chromosome chromosome --pval p_value --position base_pair_location --col y

Output: GCST90013658_buildGRCh37_plot.png 


