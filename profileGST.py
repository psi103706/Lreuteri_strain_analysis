import argparse
import os
import sys
import pandas as pd
import numpy as np

from Bio import SeqIO
from ete3 import Tree

from util import run_process

def make_GST_profile(df_sample_info, gst_profile_file_path, kraken_database_path, quiet):    
    for sample in df_sample_info.index.tolist():
        if not os.path.isfile("{0}/{1}.txt".format(gst_profile_file_path, sample)):
            proc = run_process("kraken --db {0} --threads 8 --preload --paired {1} {2} > {3}/{4}.kraken".format(kraken_database_path, df_sample_info.loc[sample, 'read_file1'], df_sample_info.loc[sample, 'read_file2'], gst_profile_file_path, sample), print_only = False)                       
            for line in iter(proc.stdout.readline, b''):
                if not quiet:
                    print(line.decode('utf-8').rstrip())
            
            proc = run_process("kraken-report --db {0} {1}/{2}.kraken > {1}/{2}.kreport".format(kraken_database_path, gst_profile_file_path, sample), print_only = False)
            for line in iter(proc.stdout.readline, b''):
                if not quiet:
                    print(line.decode('utf-8').rstrip())
                    
            proc = run_process("bracken -d {0} -i {1}/{2}.kreport -o {1}/{2}.bracken".format(kraken_database_path, gst_profile_file_path, sample), print_only = False)
            for line in iter(proc.stdout.readline, b''):
                if not quiet:
                    print(line.decode('utf-8').rstrip())
                    
def parse_bracken_result(sample, gst_profile_file_path):
    read_count = 0
    with open("{0}/{1}.kraken".format(gst_profile_file_path, sample)) as kraken_output_file:
        line = kraken_output_file.readline()
        while line != "":                       
            read_count += 1
            line = kraken_output_file.readline()
    
    df_bracken_result = pd.read_csv("{0}/{1}.bracken".format(gst_profile_file_path, sample), sep = '\t', index_col = 0)    
    df_bracken_result = df_bracken_result.loc[:, ['fraction_total_reads']]
    df_bracken_result.rename(columns = {'fraction_total_reads': 'relative abundance'}, inplace = True)
    
    return df_bracken_result, read_count

def gst_profile(df_sample_info, gst_profile_file_path):
    gst_profile = dict()
    
    for sample in df_sample_info.index.tolist():        
        if os.path.isfile("{0}/{1}.bracken".format(gst_profile_file_path, sample)) and os.path.isfile("{0}/{1}.kraken".format(gst_profile_file_path, sample)):
            df_bracken_result, read_count = parse_bracken_result(sample, gst_profile_file_path)
            if read_count >= 20000:
                gst_profile[sample] = df_bracken_result.to_dict()['relative abundance']
                    
    return pd.DataFrame(gst_profile).fillna(0.0) 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create GST profiles based on GST Kraken DB.')
    
    parser.add_argument('-i', '--input_sample_file', required = True, help = 'Input path for a sample information file.') 
    parser.add_argument('-d', '--db', required = True, help = 'GST Kraken DB.')
    parser.add_argument('-o', '--output_dir', required = True, help = 'Output directory path for GST profiles.')
    parser.add_argument('-f', '--force', help = 'Keep processing even if the output directory is already existing.', action='store_true')
    parser.add_argument('-q', '--quiet', help = 'Run quietly.', action='store_true')
    
    args = parser.parse_args()
    
    if os.path.isdir(args.output_dir):
        if not args.force:
            print('Error: Output directory already exists.')
            sys.exit()
            
    else:
        os.mkdir(args.output_dir)
    
    df_sample_info = pd.read_csv(args.input_sample_file, sep = '\t', index_col = 0)
    make_GST_profile(df_sample_info, args.output_dir, args.db, args.quiet)
    df_gst_profile = gst_profile(df_sample_info, args.output_dir)
    
    df_gst_profile.to_csv('{0}/gst_profile.txt'.format(args.output_dir), sep = '\t')