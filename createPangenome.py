import argparse
import os
import sys
import pandas as pd
import numpy as np

from Bio import SeqIO
from ete3 import Tree

from util import run_process

def generate_reuteri_pangenome(df_genome_info, pangenome_file_path, tmp_file_path, quiet):
    records = []
    for genome_id in [gn for gn in df_genome_info.index.tolist() if df_genome_info.loc[gn, 'type'] not in 'genus']:
        records = records + list(SeqIO.to_dict(SeqIO.parse("{0}".format(df_genome_info.loc[genome_id, 'cds_file']), 'fasta')).values())
        
    SeqIO.write(records, "{0}/Lactobacillus_reuteri_raw.fasta".format(pangenome_file_path), "fasta")
    
    proc = run_process("mmseqs createdb {0}/Lactobacillus_reuteri_raw.fasta {0}/Lactobacillus_reuteri_raw --dbtype 2;\n".format(pangenome_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
            
    proc = run_process("mmseqs linclust {0}/Lactobacillus_reuteri_raw {0}/Lactobacillus_reuteri_clu {1} --min-seq-id 0.9 -c 0.9 --kmer-per-seq 200 --adjust-kmer-len;\n".format(pangenome_file_path, tmp_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
            
    proc = run_process("mmseqs result2repseq {0}/Lactobacillus_reuteri_raw {0}/Lactobacillus_reuteri_clu {0}/Lactobacillus_reuteri_rep;\n".format(pangenome_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
            
    proc = run_process("mmseqs createsubdb {0}/Lactobacillus_reuteri_rep {0}/Lactobacillus_reuteri_raw_h {0}/Lactobacillus_reuteri_rep_h;\n".format(pangenome_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
            
    proc = run_process("mmseqs result2flat {0}/Lactobacillus_reuteri_raw {0}/Lactobacillus_reuteri_raw {0}/Lactobacillus_reuteri_rep {0}/Lactobacillus_reuteri.fasta;\n".format(pangenome_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
            
    proc = run_process("mmseqs createtsv {0}/Lactobacillus_reuteri_raw {0}/Lactobacillus_reuteri_raw {0}/Lactobacillus_reuteri_clu {0}/Lactobacillus_reuteri.tsv;\n".format(pangenome_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
            
def generate_lactobacillus_pangenome(df_genome_info, pangenome_file_path, tmp_file_path, quiet):         
    records = list(SeqIO.to_dict(SeqIO.parse("{0}/Lactobacillus_reuteri.fasta".format(pangenome_file_path), 'fasta')).values())
    
    for genome_id in df_genome_info[df_genome_info['type'] == 'genus'].index.tolist():
        count = 0
        for seq_id, rec in SeqIO.to_dict(SeqIO.parse("{0}".format(df_genome_info.loc[genome_id, 'cds_file']), 'fasta')).items():
            count += 1
            rec.id = rec.id + "|{0}".format(count)
            records.append(rec)
  
    SeqIO.write(records, "{0}/Lactobacillus_raw.fasta".format(pangenome_file_path), "fasta")

    proc = run_process("mmseqs createdb {0}/Lactobacillus_raw.fasta {0}/Lactobacillus_raw --dbtype 2;\n".format(pangenome_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
            
    proc = run_process("mmseqs linclust {0}/Lactobacillus_raw {0}/Lactobacillus_clu {1} --min-seq-id 0.9 -c 0.9 --kmer-per-seq 200 --adjust-kmer-len;\n".format(pangenome_file_path, tmp_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
            
    proc = run_process("mmseqs result2repseq {0}/Lactobacillus_raw {0}/Lactobacillus_clu {0}/Lactobacillus_rep;\n".format(pangenome_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
            
    proc = run_process("mmseqs createsubdb {0}/Lactobacillus_rep {0}/Lactobacillus_raw_h {0}/Lactobacillus_rep_h;\n".format(pangenome_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
    
    proc = run_process("mmseqs result2flat {0}/Lactobacillus_raw {0}/Lactobacillus_raw {0}/Lactobacillus_rep {0}/Lactobacillus.fasta;\n".format(pangenome_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
    
    proc = run_process("mmseqs createtsv {0}/Lactobacillus_raw {0}/Lactobacillus_raw {0}/Lactobacillus_clu {0}/Lactobacillus.tsv;\n".format(pangenome_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())

def parse_tsv_files(tsv_file_name, out_file_name):
    tsv_cluster_dict = dict()
    
    with open(tsv_file_name) as tsv_file:
        line = tsv_file.readline()
        while line != "":
            line = line.rstrip()
            
            info = line.split('\t')
            if info[0] not in tsv_cluster_dict.keys():
                tsv_cluster_dict[info[0]] = []
                
            tsv_cluster_dict[info[0]].append(info[1])
            
            line = tsv_file.readline()    
            
    with open(out_file_name, 'w') as out_file:
        out_file.write("Rep\tSize\tMembers\n")
        
        tsv_cluster_dict = {k: v for k, v in sorted(tsv_cluster_dict.items(), key = (lambda x: len(x[1])), reverse = True)}
        for rep, members in tsv_cluster_dict.items():
            out_file.write("{0}\t{1}\t{2}".format(rep, len(members), rep))
            
            for member in members[1:]:
                out_file.write(" {0}".format(member))
                
            out_file.write('\n')
            out_file.flush()     
            
    return tsv_cluster_dict

def get_reuteri_specific_pangenome(df_genome_info, pangenome_file_path):
    lactobacillus_cluster_dict = parse_tsv_files("{0}/Lactobacillus.tsv".format(pangenome_file_path), "{0}/Lactobacillus_cluster.txt".format(pangenome_file_path))
    
    genome_id_list = [gn for gn in df_genome_info.index.tolist() if df_genome_info.loc[gn, 'type'] not in 'genus']
    
    reuteri_specific_list = []
    reuteri_more_clustered_dict = dict()
    for rep, members in lactobacillus_cluster_dict.items():
        is_specific = True
        
        for member in members:
            info = member.split('|')
            
            genome_id = info[0]
            if "_" in genome_id:
                genome_id = genome_id.split("_")[0]
                
            if genome_id not in genome_id_list:
                is_specific = False
                break
                
        if is_specific:
            reuteri_specific_list.append(rep)
            if len(members) > 1:
                for member in members[1:]:
                    reuteri_more_clustered_dict[member] = rep
                    
    reuteri_cluster_dict = parse_tsv_files("{0}/Lactobacillus_reuteri.tsv".format(pangenome_file_path), "{0}/Lactobacillus_reuteri_cluster.txt".format(pangenome_file_path))
    
    for member, rep in reuteri_more_clustered_dict.items():
        reuteri_cluster_dict[rep] += reuteri_cluster_dict[member]
        reuteri_cluster_dict.pop(member, None)
    
    records = []
    for seq_id, record in SeqIO.to_dict(SeqIO.parse("{0}/Lactobacillus_reuteri.fasta".format(pangenome_file_path), 'fasta')).items():
        if seq_id in reuteri_cluster_dict.keys():
            records.append(record)
            
    SeqIO.write(records, "{0}/Lactobacillus_reuteri_specific.fasta".format(pangenome_file_path), "fasta")
    
    with open("{0}/Lactobacillus_reuteri_specific_cluster.txt".format(pangenome_file_path), 'w') as out_file:
        out_file.write("Rep\tSize\tMembers\n")
        
        tsv_cluster_dict = {k: v for k, v in sorted(reuteri_cluster_dict.items(), key = (lambda x: len(x[1])), reverse = True)}
        for rep, members in reuteri_cluster_dict.items():
            out_file.write("{0}\t{1}\t{2}".format(rep, len(members), rep))
            
            for member in members[1:]:
                out_file.write(" {0}".format(member))
                
            out_file.write('\n')
            out_file.flush()
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build a pan-genome database.')
    
    parser.add_argument('-i', '--input_genome_file', required = True, help = 'Input path for a genome information file.')   
    parser.add_argument('-o', '--output_dir', required = True, help = 'Output directory path for pan-genome DB.')
    parser.add_argument('-f', '--force', help = 'Keep processing even if the output directory is already existing.', action='store_true')
    parser.add_argument('-q', '--quiet', help = 'Run quietly.', action='store_true')
    
    args = parser.parse_args()
        
    if os.path.isdir(args.output_dir):
        if not args.force:
            print('Error: Output directory already exists.')
            sys.exit()
            
    else:
        os.mkdir(args.output_dir)
    
    df_genome_info = pd.read_csv(args.input_genome_file, sep = '\t', index_col = 0).fillna('')
    df_genome_info = df_genome_info[df_genome_info['cds_file'] != '']
    df_genome_info = df_genome_info.loc[[gn for gn in df_genome_info.index if df_genome_info.loc[gn, 'type'] in ['complete', 'ref', 'genus']], :]
    
    tmp_file_path = "{0}/tmp".format(args.output_dir)
    
    if not os.path.isdir(tmp_file_path):
        os.mkdir(tmp_file_path)
       
    generate_reuteri_pangenome(df_genome_info, args.output_dir, tmp_file_path, args.quiet)
    generate_lactobacillus_pangenome(df_genome_info, args.output_dir, tmp_file_path, args.quiet)
    get_reuteri_specific_pangenome(df_genome_info, args.output_dir)