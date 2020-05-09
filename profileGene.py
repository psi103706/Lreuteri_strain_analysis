import argparse
import os
import sys
import pandas as pd
import numpy as np

from Bio import SeqIO
from ete3 import Tree

from util import run_process

def prepare_reuteri_specific_db(pangenome_file_path, gene_profile_file_path, quiet):
    proc = run_process("mmseqs createdb {0}/Lactobacillus_reuteri_specific.fasta {1}/Lactobacillus_reuteri_specific --dbtype 2".format(pangenome_file_path, gene_profile_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
        
def make_gene_profile(df_sample_info, gene_profile_file_path, tmp_file_path, quiet):
    for sample in df_sample_info.index.tolist():
        if not os.path.isfile("{0}/{1}.txt".format(gene_profile_file_path, sample)):
            proc = run_process("mmseqs createdb {0} {1}/{2} --dbtype 2".format(df_sample_info.loc[sample, 'cds_file'], gene_profile_file_path, sample), print_only = False)
            for line in iter(proc.stdout.readline, b''):
                if not quiet:
                    print(line.decode('utf-8').rstrip())
            
            proc = run_process("mmseqs search {0}/{1} {0}/Lactobacillus_reuteri_specific {0}/{1}_aln {2} --min-seq-id 0.9 -c 0.9 --threads 8 --search-type 3".format(gene_profile_file_path, sample, tmp_file_path), print_only = False)
            for line in iter(proc.stdout.readline, b''):
                if not quiet:
                    print(line.decode('utf-8').rstrip())
                    
            proc = run_process("mmseqs createtsv {0}/{1} {0}/Lactobacillus_reuteri_specific {0}/{1}_aln {0}/{1}.txt --threads 8".format(gene_profile_file_path, sample), print_only = False)
            for line in iter(proc.stdout.readline, b''):
                if not quiet:
                    print(line.decode('utf-8').rstrip())
                    
def parse_mmseqs_result(sample, gene_profile_file_path):
    hit_dict = dict()
    
    with open('{0}/{1}.txt'.format(gene_profile_file_path, sample)) as mmseqs_result_file:
        line = mmseqs_result_file.readline()
        while line != "":
            line = line.rstrip()            
            info = line.split('\t')
            
            query = info[0]
            target = info[1]
            identity = float(info[3])
            
            if query in hit_dict.keys():
                if identity > hit_dict[query][1]:
                    hit_dict[query] = [target, identity]
                    
            else:
                hit_dict[query] = [target, identity]
            
            line = mmseqs_result_file.readline()
    
    total_count = 0
    
    gene_profile_dict = dict()
    for query, info in hit_dict.items():
        total_count += 1
        
        try:
            gene_profile_dict[info[0]] += 1
            
        except:
            gene_profile_dict[info[0]] = 1            
            
    return gene_profile_dict, hit_dict, total_count

def gene_profile(df_sample_info, gene_profile_file_path):
    gene_profile = dict()    

    for sample in df_sample_info.index.tolist():
        if os.path.isfile("{0}/{1}.txt".format(gene_profile_file_path, sample)):
            gene_profile_dict, _, total_count = parse_mmseqs_result(sample, gene_profile_file_path)
            if total_count >= 500:
                gene_profile[sample] = gene_profile_dict
                    
    return pd.DataFrame(gene_profile).fillna(0.0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create gene profiles based on Pan-genome DB.')
    
    parser.add_argument('-i', '--input_sample_file', required = True, help = 'Input path for a sample information file.') 
    parser.add_argument('-d', '--db', required = True, help = 'Pan-genome DB.')
    parser.add_argument('-o', '--output_dir', required = True, help = 'Output directory path for gene profiles.')
    parser.add_argument('-f', '--force', help = 'Keep processing even if the output directory is already existing.', action='store_true')
    parser.add_argument('-q', '--quiet', help = 'Run quietly.', action='store_true')
    
    args = parser.parse_args()
    
    if os.path.isdir(args.output_dir):
        if not args.force:
            print('Error: Output directory already exists.')
            sys.exit()
            
    else:
        os.mkdir(args.output_dir)
    
    tmp_file_path = "{0}/tmp".format(args.output_dir)
    
    if not os.path.isdir(tmp_file_path):
        os.mkdir(tmp_file_path)
    
    df_sample_info = pd.read_csv(args.input_sample_file, sep = '\t', index_col = 0)
    prepare_reuteri_specific_db(args.db, args.output_dir, args.quiet)
    make_gene_profile(df_sample_info, args.output_dir, tmp_file_path, args.quiet)
    
    df_gene_profile = gene_profile(df_sample_info, args.output_dir)
    
    df_gene_profile.to_csv('{0}/gene_profile.txt'.format(args.output_dir), sep = '\t')