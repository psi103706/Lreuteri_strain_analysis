import argparse
import os
import sys
import pandas as pd
import numpy as np

from Bio import SeqIO
from ete3 import Tree

from util import run_process

def extract_core_gene(df_genome_info, output_dir, core_gene_file_path, quiet):      
    df_complete = df_genome_info[df_genome_info['type'] == 'complete']  
    
    for i in range(df_complete.index.size):
        genome_id = df_complete.index.tolist()[i]
        proc = run_process("prokka {0} --outdir {1}/gff --prefix {2} --force".format(df_complete.loc[genome_id, 'genome_file'], output_dir, i + 1), print_only = False)
        for line in iter(proc.stdout.readline, b''):
            if not quiet:
                print(line.decode('utf-8').rstrip())           
    
    proc = run_process("roary {0}/gff/*.gff -f {0}/roary_output -p 8 -e".format(output_dir), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
     
    core_gene_name_list = write_core_gene_file("{0}/roary_output".format(output_dir), core_gene_file_path, df_complete)
    return core_gene_name_list

def write_core_gene_file(roary_output_path, core_gene_file_path, df_complete):
    df_gene_pre_abs = pd.read_csv("{0}/gene_presence_absence.csv".format(roary_output_path), index_col = 0)
    core_gene_list = df_gene_pre_abs[(df_gene_pre_abs['No. isolates'] == df_complete.index.size) & (df_gene_pre_abs['No. sequences'] == df_complete.index.size) & (df_gene_pre_abs['Avg sequences per isolate'] == 1)].index.tolist()
    
    core_gene_name_list = []
    core_gene_seq_list = []
    for record in SeqIO.parse("{0}/pan_genome_reference.fa".format(roary_output_path), "fasta") :
        if record.description.split(' ')[1] in core_gene_list and len(record.seq) > 100:
            record.id = record.description.split(' ')[1]
            record.description = ''
            
            core_gene_name_list.append(record.id)
            core_gene_seq_list.append(record)
    
    if not os.path.isfile("{0}/core_gene_ref.fasta".format(core_gene_file_path)):
        SeqIO.write(core_gene_seq_list, "{0}/core_gene_ref.fasta".format(core_gene_file_path), "fasta")
    
    return sorted(core_gene_name_list)

def search_core_gene(df_genome_info, genome_id, core_gene_file_path, delta_output_path, snp_output_path, quiet):
    if not os.path.isfile("{0}/{1}.delta".format(delta_output_path, genome_id)):
        print("{0}/{1}.delta".format(delta_output_path, genome_id))
        proc = run_process("nucmer {0}/core_gene_ref.fasta {1} -p {2}/{3}".format(core_gene_file_path, df_genome_info.loc[genome_id, 'genome_file'], delta_output_path, genome_id), print_only = False)

        for line in iter(proc.stdout.readline, b''):
            if not quiet:
                print(line.decode('utf-8').rstrip())   
    
    if not os.path.isfile("{0}/{1}_filter.delta".format(delta_output_path, genome_id)):
        proc = run_process("delta-filter -r -q {0}/{1}.delta > {0}/{1}_filter.delta".format(delta_output_path, genome_id), print_only = False)

        for line in iter(proc.stdout.readline, b''):
            if not quiet:
                print(line.decode('utf-8').rstrip())

    core_block_dict = dict()
    core_len_dict = dict()        
    
    with open("{0}/{1}_filter.delta".format(delta_output_path, genome_id)) as delta_file:
        line = delta_file.readline()        
        while line != "":
            line = line.rstrip()

            if line[0] == '>':
                info = line[1:].split(' ')
                core, ref_len = info[0], int(info[2])
                core_len_dict[core] = ref_len

                if core not in core_block_dict.keys():
                    core_block_dict[core] = []

            elif len(line.split(' ')) == 7:
                core_block_dict[core].append((int(line.split(' ')[0]), int(line.split(' ')[1])))            
                
            line = delta_file.readline()
                
    core_gene_list = sorted(core_len_dict.keys())    
   
    core_empty_block_dict = dict()    
    for core in core_block_dict.keys():              
        block_list = sorted(core_block_dict[core], key=lambda block: block[0])
        length = core_len_dict[core]
        
        curr_start = 1
        for start, end in block_list:                
            if start - 1 >= curr_start:
                if core not in core_empty_block_dict.keys():
                    core_empty_block_dict[core] = []
                    
                core_empty_block_dict[core].append((curr_start, start - 1))
            
            if curr_start <= end:
                curr_start = end + 1
        
        if length >= curr_start:
            if core not in core_empty_block_dict.keys():
                core_empty_block_dict[core] = []
                
            core_empty_block_dict[core].append((curr_start, length))            

   
    core_gap_dict = dict()
    core_var_dict = dict()
    if not os.path.isfile("{0}/{1}.snp".format(snp_output_path, genome_id)):
        with open("{0}/{1}.snp".format(snp_output_path, genome_id), 'w') as snp_file:
            proc = run_process("show-snps -H -T {0}/{1}_filter.delta".format(delta_output_path, genome_id), print_only = False)

            for line in iter(proc.stdout.readline, b''):       
                line = line.decode('utf-8').rstrip()
                if not quiet:
                    print(line)
                snp_file.write(line + '\n')
                
    with open("{0}/{1}.snp".format(snp_output_path, genome_id)) as snp_file:
        line = snp_file.readline()            
        while line != "":
            line = line.rstrip()                
            info = line.split('\t')
                
            pos = int(info[0])
            ref_char = '-' if info[1] == '.' else info[1]
            query_char = '-' if info[2] == '.' else info[2]
            core = info[10]
                
            if ref_char == '-':
                if pos in core_gap_dict.keys():
                    if core not in core_gap_dict.keys():
                        core_gap_dict[core] = dict()
                
                    core_gap_dict[core][pos].append(query_char)
                        
                else:
                    if core not in core_gap_dict.keys():
                        core_gap_dict[core] = dict()
                        
                    core_gap_dict[core][pos] = [query_char]
                        
            else:
                if core not in core_var_dict.keys():
                    core_var_dict[core] = dict()
                    
                core_var_dict[core][pos] = (ref_char, query_char)
                
            line = snp_file.readline()  
    
    return core_gene_list, core_empty_block_dict, core_gap_dict, core_var_dict

def parse_align_fasta(df_genome_info, genome_id, core_gene_file_path, delta_output_path, snp_output_path, quiet):    
    core_gene_list, core_empty_block_dict, core_gap_dict, core_var_dict = search_core_gene(df_genome_info, genome_id, core_gene_file_path, delta_output_path, snp_output_path, quiet)
    
    if not os.path.isfile("{0}/{1}_aln.fasta".format(core_gene_file_path, genome_id)):     
        with open("{0}/{1}_aln.fasta".format(core_gene_file_path, genome_id), 'w') as aln_file:
            ref_seq_record_dict = SeqIO.to_dict(SeqIO.parse("{0}/core_gene_ref.fasta".format(core_gene_file_path), 'fasta'))
            
            for core in core_gene_list:
                empty_block_list = core_empty_block_dict[core] if core in core_empty_block_dict.keys() else None
                gap_dict = core_gap_dict[core] if core in core_gap_dict.keys() else None
                var_dict = core_var_dict[core] if core in core_var_dict.keys() else None
                
                query_seq_align = ""
                ref_seq_align = ""
                
                if empty_block_list is not None or gap_dict is not None or var_dict is not None:
                    empty_block_dict = None
                    if empty_block_list is not None:
                        empty_block_dict = dict()
                        
                        for i in range(len(empty_block_list)):
                            empty_block_dict[empty_block_list[i][0]] = i
                    
                    if gap_dict is not None and 0 in gap_dict.keys():
                        q_seq = ''.join(gap_dict[0])
                        
                        query_seq_align = q_seq
                        ref_seq_align = '-' * len(q_seq)
                    
                    r_seq = str(ref_seq_record_dict[core].seq)
                    
                    curr = 1
                    while curr - 1 < len(r_seq):
                        if empty_block_dict is not None and curr in empty_block_dict.keys():
                            empty_block = empty_block_list[empty_block_dict[curr]]
                            
                            query_seq_align += '-' * (empty_block[1] - empty_block[0] + 1)
                            ref_seq_align += r_seq[empty_block[0] - 1 : empty_block[1]]
                            
                            curr = empty_block[1] + 1
                        
                        else:
                            if var_dict is not None and curr in var_dict.keys():
                                query_seq_align += var_dict[curr][1]
                                ref_seq_align += var_dict[curr][0]
                                
                                if r_seq[curr - 1] != var_dict[curr][0]:
                                    print(r_seq, curr, r_seq[curr - 1], var_dict[curr][0], var_dict[curr][1])
                                
                            else:
                                query_seq_align += r_seq[curr - 1]
                                ref_seq_align += r_seq[curr - 1]
                            
                            if gap_dict is not None and curr in gap_dict.keys():
                                q_seq = ''.join(gap_dict[curr])

                                query_seq_align += q_seq
                                ref_seq_align += '-' * len(q_seq)
                        
                            curr += 1
                            
                else:
                    query_seq_align = str(ref_seq_record_dict[core].seq)
                    ref_seq_align = str(ref_seq_record_dict[core].seq)
                    
                aln_file.write(">{0}_q\n{1}\n>{0}_r\n{2}\n".format(core, query_seq_align, ref_seq_align))
    
    core_align_query_ref_dict = dict()
    
    aln_seq_record_dict = SeqIO.to_dict(SeqIO.parse("{0}/{1}_aln.fasta".format(core_gene_file_path, genome_id), 'fasta'))
    for core in core_gene_list:        
        core_align_query_ref_dict[core] = (str(aln_seq_record_dict[core + "_q"].seq), str(aln_seq_record_dict[core + "_r"].seq))
        
    return core_align_query_ref_dict

def calc_total_gap_count(df_genome_info, ref_id, core_gene_name_list, core_gene_file_path, delta_output_path, snp_output_path, quiet):    
    gap_count_dict = dict()
    
    for i in range(df_genome_info.index.size):
        genome_id = df_genome_info.index.tolist()[i]
        
        core_align_query_ref_dict = parse_align_fasta(df_genome_info, genome_id, core_gene_file_path, delta_output_path, snp_output_path, quiet)
            
        for core, seq_pair in core_align_query_ref_dict.items():
            if core not in gap_count_dict.keys():
                gap_count_dict[core] = dict()
                
            nuc_count = 0
            gap_count = 0
            for pos in range(len(seq_pair[1])):
                if seq_pair[1][pos] != '-':
                    if gap_count > 0:
                        if nuc_count not in gap_count_dict[core].keys():
                            gap_count_dict[core][nuc_count] = (gap_count, genome_id)
                                
                        elif gap_count_dict[core][nuc_count][0] < gap_count:
                            gap_count_dict[core][nuc_count] = (gap_count, genome_id)

                        gap_count = 0
                            
                    nuc_count += 1
                        
                else:
                    gap_count += 1
                    
            if gap_count > 0:
                if nuc_count not in gap_count_dict[core].keys():
                    gap_count_dict[core][nuc_count] = (gap_count, genome_id)
                                
                elif gap_count_dict[core][nuc_count][0] < gap_count:
                    gap_count_dict[core][nuc_count] = (gap_count, genome_id)
                     
    return gap_count_dict 

def merge_gap_for_msa(df_genome_info, ref_id, core_gene_name_list, core_gene_file_path, delta_output_path, snp_output_path, quiet):   
    gap_count_dict = calc_total_gap_count(df_genome_info, ref_id, core_gene_name_list, core_gene_file_path, delta_output_path, snp_output_path, quiet)
    
    partition_dict = dict()
        
    with open("{0}/all_align_core.fasta".format(core_gene_file_path), 'w') as all_align_core_file:
        core_align_query_ref_dict = parse_align_fasta(df_genome_info, ref_id, core_gene_file_path, delta_output_path, snp_output_path, quiet)  
        core_len_dict = dict()
                
        all_align_core_file.write(">{0}\n".format(ref_id))
        
        start_pos = 1
        for core in core_gene_name_list:
            ref_seq = core_align_query_ref_dict[core][1]            
            ref_seq_add = ""            
            
            if len(gap_count_dict[core]) == 0:
                ref_seq_add = ref_seq
            
            else:
                nuc_count = 0
                gap_checked = False
                
                for i in range(len(ref_seq)): 
                    if nuc_count in gap_count_dict[core]:
                        if not gap_checked:
                            for j in range(gap_count_dict[core][nuc_count][0]):
                                ref_seq_add += '-'
                                
                            gap_checked = True
                    
                    if ref_seq[i] != '-':                                
                        ref_seq_add += ref_seq[i]
                        nuc_count += 1
                        
                        if gap_checked:
                            gap_checked = False
                            
                if nuc_count in gap_count_dict[core]:
                    if not gap_checked:
                        for j in range(gap_count_dict[core][nuc_count][0]):
                            ref_seq_add += '-'
                                
                        gap_checked = True

            all_align_core_file.write(ref_seq_add)                    
            core_len_dict[core] = len(ref_seq_add)
            
            partition_dict[core] = [start_pos, start_pos + len(ref_seq_add) - 1]
            start_pos += len(ref_seq_add)
            
        all_align_core_file.write("\n")

        for i in range(df_genome_info.index.size):
            genome_id = df_genome_info.index.tolist()[i]

            if genome_id != ref_id:                
                core_align_query_ref_dict = parse_align_fasta(df_genome_info, genome_id, core_gene_file_path, delta_output_path, snp_output_path, quiet)
                
                all_align_core_file.write(">{0}\n".format(genome_id))                
                
                for core in core_gene_name_list:
                    query_seq_add = ""
                    
                    if core in core_align_query_ref_dict.keys():
                        ref_seq = core_align_query_ref_dict[core][1]      
                        query_seq = core_align_query_ref_dict[core][0]                                    
                        
                        if len(gap_count_dict[core]) == 0:
                            query_seq_add = query_seq
                            
                        else:                          
                            nuc_count = 0
                            gap_checked = False
                            
                            for i in range(len(ref_seq)):
                                if nuc_count in gap_count_dict[core]:
                                    if not gap_checked:                            
                                        gap_count_total = gap_count_dict[core][nuc_count][0]
                                        gap_count = 0
                                    
                                        j = i
                                        while j < len(ref_seq) and ref_seq[j] == "-":
                                            gap_count += 1
                                            j += 1
                                        
                                        query_seq_add += query_seq[i : i + gap_count]
                                        for k in range(gap_count_total - gap_count):
                                            query_seq_add += "-"
                                            
                                        gap_checked = True
                                            

                                if ref_seq[i] != '-':                                                                       
                                    query_seq_add += query_seq[i]
                                    nuc_count += 1
                                    
                                    if gap_checked:
                                        gap_checked = False
                            
                            if nuc_count in gap_count_dict[core]:
                                if not gap_checked:                            
                                    gap_count_total = gap_count_dict[core][nuc_count][0]
                                    gap_count = 0
                                    
                                    j = i
                                    while j < len(ref_seq) and ref_seq[j] == "-":
                                        gap_count += 1
                                        j += 1
                                        
                                    query_seq_add += query_seq[i : i + gap_count]
                                    for k in range(gap_count_total - gap_count):
                                        query_seq_add += "-"
                                            
                                    gap_checked = True
                                        
                        for i in range(core_len_dict[core] - len(query_seq_add)):
                            query_seq_add += "-"
                        
                    else:
                        for i in range(core_len_dict[core]):
                            query_seq_add += "-"
                    
                    if len(query_seq_add) != core_len_dict[core]:
                        print(genome_id, core, (core in core_align_query_ref_dict.keys()), len(query_seq_add), core_len_dict[core])
                        
                    all_align_core_file.write(query_seq_add)
                    
                all_align_core_file.write("\n")
                
        with open("{0}/partition.txt".format(core_gene_file_path), 'w') as partition_file:
            for gene, pos in partition_dict.items():
                partition_file.write("DNA, {0} = {1}-{2}\n".format(gene, pos[0], pos[1]))

def generate_ml_tree(suffix, core_gene_file_path, quiet):
    if not os.path.isfile("{1}/reuteri_core{0}.nwk".format("" if suffix is None else ("_" + suffix), core_gene_file_path)):
        proc = run_process("raxmlHPC-PTHREADS -s {0}/all_align_core.fasta -n reuteri_core{1} -T 8 -m GTRCAT -p 123".format(core_gene_file_path, "" if suffix is None else ("_" + suffix)), print_only = False)

        for line in iter(proc.stdout.readline, b''):
            if line.strip() and not quiet:
                print(line.decode('utf-8').rstrip())

        run_process("mv RAxML_bestTree.reuteri_core{0} {1}/reuteri_core{0}.nwk".format("" if suffix is None else ("_" + suffix), core_gene_file_path))
        run_process("rm RAxML_*.reuteri_core{0}*".format("" if suffix is None else ("_" + suffix)))
    
    return Tree("{1}/reuteri_core{0}.nwk".format("" if suffix is None else ("_" + suffix), core_gene_file_path))

def get_final_tree(tree_all, df_genome_info, core_gene_file_path, out_tree_file_name):
    if not os.path.isfile("{0}/{1}".format(core_gene_file_path, out_tree_file_name)):
        outgroup_node = None
        outgroup_list = [gn for gn in df_genome_info.index.tolist() if df_genome_info.loc[gn, 'type'] == 'out']
        if len(outgroup_list) == 1:
            outgroup_node = tree_all&outgroup_list[0]            
            
        else:
            outgroup_node = tree_all.get_common_ancestor(outgroup_list)
            
        tree_all.set_outgroup(outgroup_node)

        outgroup_node.detach()
        tree_all.get_tree_root().get_children()[0].delete()
        tree_all.ladderize()

        tree_all.write(outfile = "{0}/{1}".format(core_gene_file_path, out_tree_file_name))
        return tree_all
    
    else:
        return Tree("{0}/{1}".format(core_gene_file_path, out_tree_file_name), format = 1)

def calc_snp_diff(df_genome_info, genome_id1, genome_id2, core_gene_file_path, delta_output_path, snp_output_path, quiet):
    core_align_query_ref_dict1 = parse_align_fasta(df_genome_info, genome_id1, core_gene_file_path, delta_output_path, snp_output_path, quiet)
    core_align_query_ref_dict2 = parse_align_fasta(df_genome_info, genome_id2, core_gene_file_path, delta_output_path, snp_output_path, quiet)
    
    diff_cnt = 0
    diff_no_gap_count = 0
    
    for core in core_align_query_ref_dict1.keys():
        if core in core_align_query_ref_dict2.keys():      
            query1 = core_align_query_ref_dict1[core][0]
            query2 = core_align_query_ref_dict2[core][0]

            ref1 = core_align_query_ref_dict1[core][1]
            ref2 = core_align_query_ref_dict2[core][1]
            
            new_q1 = ""
            new_q2 = ""
            
            gap_map1 = dict(); gap_map2 = dict()
            curr_pos1 = 0; curr_pos2 = 0;
            
            while curr_pos1 < len(ref1) and curr_pos2 < len(ref2):
                if ref1[curr_pos1] == ref2[curr_pos2]:                   
                    new_q1 = new_q1 + query1[curr_pos1]
                    new_q2 = new_q2 + query2[curr_pos2]
                    
                    curr_pos1 += 1
                    curr_pos2 += 1
                    
                elif ref1[curr_pos1] != '-' and ref2[curr_pos2] == '-':
                    new_q1 = new_q1 + '-'
                    new_q2 = new_q2 + query2[curr_pos2]
                    
                    curr_pos2 += 1
                    
                elif ref1[curr_pos1] == '-' and ref2[curr_pos2] != '-':
                    new_q1 = new_q1 + query1[curr_pos1]
                    new_q2 = new_q2 + '-'
                    
                    curr_pos1 += 1
                    
            if curr_pos1 < len(ref1) and curr_pos2 == len(ref2):
                new_q2 = new_q2 + ('-' * (len(ref1) - curr_pos1))
                
            
            if curr_pos1 == len(ref1) and curr_pos2 < len(ref2):
                new_q1 = new_q1 + ('-' * (len(ref2) - curr_pos2))
                
            for c1, c2 in zip(new_q1, new_q2):
                if c1 != c2:
                    diff_cnt += 1
                    
                    if c1 != '-' and c2 != '-':
                        diff_no_gap_count += 1
                        
    return diff_cnt, diff_no_gap_count

    
def get_snp_matrix(df_genome_info, core_gene_file_path, delta_output_path, snp_output_path, quiet):
    uid_list = [gn for gn in df_genome_info.index.tolist() if df_genome_info.loc[gn, 'type'] not in 'out']
    
    snp_matrix = np.zeros((len(uid_list), len(uid_list)))
    snp_matrix_no_gap = np.zeros((len(uid_list), len(uid_list)))
    
    if not os.path.isfile("{0}/snp_matrix.txt".format(core_gene_file_path)):
        with open("{0}/snp_matrix.txt".format(core_gene_file_path), 'w') as snp_matrix_file:
            with open("{0}/snp_matrix_no_gap.txt".format(core_gene_file_path), 'w') as snp_matrix_no_gap_file:
                for uid in uid_list:
                    snp_matrix_file.write("\t{0}".format(uid))
                    snp_matrix_no_gap_file.write("\t{0}".format(uid))
                    
                snp_matrix_file.write("\n")    
                snp_matrix_no_gap_file.write("\n")   

                for i in range(len(uid_list)):
                    uid_i = str(uid_list[i])
                    snp_matrix_file.write(uid_i)
                    snp_matrix_no_gap_file.write(uid_i)

                    for j in range(len(uid_list)):
                        uid_j = str(uid_list[j])           

                        count = 0
                        count_no_gap = 0
                        if i == j:
                            pass

                        elif i < j:
                            diff_cnt, diff_no_gap_count = calc_snp_diff(df_genome_info, uid_i, uid_j, core_gene_file_path, delta_output_path, snp_output_path, quiet)
                            
                            snp_matrix[i][j] = diff_cnt
                            snp_matrix_no_gap[i][j] = diff_no_gap_count

                            count = diff_cnt
                            count_no_gap = diff_no_gap_count

                        else:
                            count = snp_matrix[j][i]
                            count_no_gap = snp_matrix_no_gap[j][i]

                        snp_matrix_file.write("\t{0}".format(count))
                        snp_matrix_no_gap_file.write("\t{0}".format(count_no_gap))
                        
                    snp_matrix_file.write("\n")
                    snp_matrix_file.flush()
                    
                    snp_matrix_no_gap_file.write("\n")
                    snp_matrix_no_gap_file.flush()
            
    df_snp = pd.read_csv("{0}/snp_matrix.txt".format(core_gene_file_path), sep = '\t', index_col = 0)
    df_snp.index = df_snp.index.map(str)
    df_snp.columns = df_snp.columns.map(str)
    
    return df_snp 

def cluster_strain_types(df_genome_info, tree, threshold, core_gene_file_path, delta_output_path, snp_output_path, quiet):
    uid_list = [nd.name for nd in tree.get_leaves()]
    df_snp = get_snp_matrix(df_genome_info, core_gene_file_path, delta_output_path, snp_output_path, quiet)
    
    sorted_uid_list = sorted(uid_list)
    used_uid_set = set()
    
    cluster_index = 0
    cluster_map = dict()
    
    for uid in sorted_uid_list:
        if uid not in used_uid_set:
            cluster_index += 1
            
            used_uid_set.add(uid)
            cluster_map[cluster_index] = [uid]            
            
            node = tree.search_nodes(name = uid)[0]
            while not node.is_root():
                node = node.up
                leaf_uid_list = [leaf_node.name for leaf_node in node.get_leaves()]
                
                if df_snp.loc[leaf_uid_list, leaf_uid_list].max().max() < threshold:
                    cluster_map[cluster_index] = leaf_uid_list
                    
                    for leaf_uid in leaf_uid_list:
                        used_uid_set.add(leaf_uid)                    
                
                else:
                    break         
    
    sorted_cluster_map = {i + 1: members for i, members in enumerate(sorted(cluster_map.values(), key = lambda x: uid_list.index(x[0])))}
    
    subcluster_map = dict()
    for cluster_index, members in sorted_cluster_map.items():
        subcluster = []
        
        sorted_uid_list_ = sorted(members)
        used_uid_set_ = set()
        
        for uid in sorted_uid_list_:
            if uid not in used_uid_set_:                
                subcluster_ = set()
                
                used_uid_set_.add(uid)    
                subcluster_.add(uid)
                
                node = tree.search_nodes(name = uid)[0]
                while node != tree.get_common_ancestor(sorted_uid_list_):
                    node = node.up
                    leaf_uid_list = [leaf_node.name for leaf_node in node.get_leaves()]

                    if df_snp.loc[leaf_uid_list, leaf_uid_list].max().max() < (threshold / 2):

                        for leaf_uid in leaf_uid_list:
                            used_uid_set_.add(leaf_uid)    
                            subcluster_.add(leaf_uid)

                    else:
                        break
                        
                subcluster_ = sorted(list(subcluster_))
                subcluster.append(subcluster_)            
         
        subcluster_map[cluster_index] = sorted(subcluster, key = lambda x: uid_list.index(x[0]))
        
    uid_type_dict = {member: "GST_{0}{1}".format(group_id, "" if len(subgroups) == 1 else "_{0}".format(i + 1)) for group_id, subgroups in subcluster_map.items() for i in range(len(subgroups)) for member in subgroups[i]}
    df_snp.loc[uid_list, uid_list].to_csv("{0}/snp_matrix.txt".format(core_gene_file_path), sep = '\t')
    
    return subcluster_map, uid_type_dict, df_snp.loc[uid_list, uid_list]

def annotate_internal_tree(tree_file_name, subgroup_dict, core_gene_file_path):
    tree = Tree("{0}/{1}".format(core_gene_file_path, tree_file_name), format = 0)    
    
    subgroup_name_dict = dict()
    for index, subgroups in subgroup_dict.items():
        subgroup_name_dict["GST_{0}".format(index)] = []
        
        if len(subgroups) > 1:
            for subgroup_ind in range(len(subgroups)):
                subgroup = subgroups[subgroup_ind]

                subgroup_nd = None
                if len(subgroup) > 1:
                    subgroup_nd = tree.get_common_ancestor(subgroup)

                else:
                    subgroup_nd = tree&subgroup[0]

                subgroup_nd.name = "GST_{0}_{1}".format(index, subgroup_ind + 1)
                subgroup_name_dict["GST_{0}".format(index)].append(subgroup_nd.name)

            group_nd = tree.get_common_ancestor(subgroup_name_dict["GST_{0}".format(index)])
            group_nd.name = "GST_{0}".format(index)
            
        else:
            group_nd = None
            
            if len(subgroups[0]) > 1:
                group_nd = tree.get_common_ancestor(subgroups[0])

            else:
                group_nd = tree&subgroups[0][0]
                
            group_nd.name = "GST_{0}".format(index)
            
            subgroup_name_dict["GST_{0}".format(index)].append(group_nd.name)
    
    tree.prune([gst_sub_ for gst, gst_sub in subgroup_name_dict.items() for gst_sub_ in gst_sub])
    
    for nd in tree.traverse():
        if nd.name == '':
            gst_set = set()
            for nd_leaf in nd.get_leaves():
                gst_set.add(nd_leaf.name.split("_")[1])
                
            if len(gst_set) == 1:
                gst_sub_list = []
                for nd_leaf in nd.get_leaves():
                    gst_sub_list.append(int(nd_leaf.name.split("_")[2]))
                
                nd.name = "GST_{0}_{1}_{2}".format(list(gst_set)[0], min(gst_sub_list), max(gst_sub_list))
            
            else:
                nd.name = "INTERNAL_{0}_{1}".format(min(list(gst_set)), max(list(gst_set)))
    
    
    tree.write(outfile = "{0}/reuteri_core_rep_internal.nwk".format(core_gene_file_path), format = 1)

def generate_kraken_taxonomy_file(tree_rep_interval, kraken_taxonomy_file_path):
    tail_node_dmp = "	|		|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	1	|		|"
    tail_name_dmp = "	|		|	scientific name	|"
    
    interval_str = "	|	"
    
    with open("{0}/nodes.dmp".format(kraken_taxonomy_file_path), 'w') as nodes_file:
        with open("{0}/names.dmp".format(kraken_taxonomy_file_path), 'w') as names_file:
            nodes_file.write("1{0}1{0}no rank{1}\n".format(interval_str, tail_node_dmp))
            nodes_file.write("812947{0}1{0}family{1}\n".format(interval_str, tail_node_dmp))
            
            names_file.write("1{0}root{1}\n".format(interval_str, tail_name_dmp))
            names_file.write("812947{0}Lactobacillus reuteri{1}\n".format(interval_str, tail_name_dmp))
           
            internal_index = 1
            
            name_id_dict = {'': '812947'}
            type_internal_index_dict = dict()
            
            for nd in tree_rep_interval.traverse(strategy = 'preorder'):
                if nd.name == '':
                    continue
                
                else:
                    id_ = None
                    
                    if nd.name.startswith("INTERNAL"):
                        id_ = "10000{0:02}".format(internal_index)
                        internal_index += 1
                        
                        nodes_file.write("{0}{1}{2}{1}genus{3}\n".format(id_, interval_str, name_id_dict[nd.up.name], tail_node_dmp))
                        names_file.write("{0}{1}{2}{3}\n".format(id_, interval_str, nd.name, tail_name_dmp))
                    
                    elif nd.name.startswith("GST"):
                        name_info = nd.name.split('_')
                        type_id = int(name_info[1])
                        
                        if len(name_info) == 2:
                            id_ = "9000{0:02}".format(type_id)
                            nodes_file.write("{0}{1}{2}{1}species{3}\n".format(id_, interval_str, name_id_dict[nd.up.name], tail_node_dmp))
                            names_file.write("{0}{1}{2}{3}\n".format(id_, interval_str, nd.name, tail_name_dmp))
                            
                            type_internal_index_dict[type_id] = 1
                            
                        elif len(name_info) == 3:
                            sub_type_id = int(name_info[2])
                            id_ = "9{0:02}0{1:02}".format(sub_type_id, type_id)
                            nodes_file.write("{0}{1}{2}{1}strain{3}\n".format(id_, interval_str, name_id_dict[nd.up.name], tail_node_dmp))
                            names_file.write("{0}{1}{2}{3}\n".format(id_, interval_str, nd.name, tail_name_dmp))
                            
                        elif len(name_info) == 4:  
                            id_ = "11{0:02}0{1:02}".format(type_internal_index_dict[type_id], type_id)
                            nodes_file.write("{0}{1}{2}{1}strain{3}\n".format(id_, interval_str, name_id_dict[nd.up.name], tail_node_dmp))
                            names_file.write("{0}{1}{2}{3}\n".format(id_, interval_str, nd.name, tail_name_dmp))
                            
                            type_internal_index_dict[type_id] += 1
                            
                    name_id_dict[nd.name] = id_
                    
    return name_id_dict

def generate_kraken_fasta_file(subgroup_dict, name_id_dict, core_gene_file_path, kraken_file_path, kraken_fasta_file_path):
    kraken_fasta_tag = "|kraken:taxid|"
    
    seq_to_id_dict = dict()
    for index, subgroups in subgroup_dict.items():
        if len(subgroups) > 1:
            for i in range(len(subgroups)):
                subgroup_rep_uid = subgroups[i][0]          
                node_name = "GST_{0}_{1}".format(index, i + 1)
                
                records = []
                for seq_id, record in SeqIO.to_dict(SeqIO.parse("{0}/{1}_aln.fasta".format(core_gene_file_path, subgroup_rep_uid), 'fasta')).items():
                    if record.id.endswith("_q"):
                        record.id = record.id[:-2]
                        record.id = record.id + kraken_fasta_tag + name_id_dict[node_name]
                        record.description = ""

                        records.append(record)
                        seq_to_id_dict[record.id] = name_id_dict[node_name]
                
                SeqIO.write(records, "{0}/{1}.fasta".format(kraken_fasta_file_path, subgroup_rep_uid), "fasta")
                
        else:
            group_rep_uid = subgroups[0][0]          
            node_name = "GST_{0}".format(index)
            
            records = []
            for seq_id, record in SeqIO.to_dict(SeqIO.parse("{0}/{1}_aln.fasta".format(core_gene_file_path, group_rep_uid), 'fasta')).items():
                if record.id.endswith("_q"):
                    record.id = record.id[:-2]
                    record.id = record.id + kraken_fasta_tag + name_id_dict[node_name]
                    record.description = ""

                    records.append(record)
                    seq_to_id_dict[record.id] = name_id_dict[node_name]
            
            SeqIO.write(records, "{0}/{1}.fasta".format(kraken_fasta_file_path, group_rep_uid), "fasta")
            
    with open("{0}/seqid2taxid.map".format(kraken_file_path), 'w') as seq_to_tax_file:
        for seq_id, tax_id in seq_to_id_dict.items():
            seq_to_tax_file.write("{0}\t{1}\n".format(seq_id, tax_id))    

def build_kraken_bracken_db(subgroup_dict, kraken_file_path, kraken_fasta_file_path, quiet):
    for index, subgroups in subgroup_dict.items():        
        for subgroup in subgroups:
            proc = run_process("kraken-build --db {0} --add-to-library {1}/{2}.fasta;\n".format(kraken_file_path, kraken_fasta_file_path, subgroup[0]), print_only = False)
            for line in iter(proc.stdout.readline, b''):
                if not quiet:
                    print(line.decode('utf-8').rstrip())
                
    proc = run_process("kraken-build --db {0} --build --threads 8;\n".format(kraken_file_path), print_only = False)
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
            
    proc = run_process("bracken-build -k 31 -d {0} -x /jc/bin/kraken/ -t 8;\n".format(kraken_file_path), print_only = False) 
    for line in iter(proc.stdout.readline, b''):
        if not quiet:
            print(line.decode('utf-8').rstrip())
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Pick GSTs from complete reference genomes and build a GST Kraken database.')
    
    parser.add_argument('-i', '--input_genome_file', required = True, help = 'Input path for a genome information file.')
    parser.add_argument('-o', '--output_dir', required = True, help = 'Output directory path for GST picking.')
    parser.add_argument('-d','--db', required = True, help = 'Output directory path for GST Kraken DB.')
    parser.add_argument('--snv_threshold', type = int, default = 20000, help = 'The number of SNVs used to cluster GST.')
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
    df_genome_info = df_genome_info[df_genome_info['genome_file'] != '']
    df_genome_info = df_genome_info.loc[[gn for gn in df_genome_info.index if df_genome_info.loc[gn, 'type'] in ['complete', 'ref', 'out']], :]
    
    core_gene_file_path = "{0}/core_gene".format(args.output_dir)
    delta_output_path = "{0}/delta".format(args.output_dir)
    snp_output_path = "{0}/snp".format(args.output_dir)
    
    if not os.path.isdir(core_gene_file_path):
        os.mkdir(core_gene_file_path)
        
    if not os.path.isdir(delta_output_path):
        os.mkdir(delta_output_path)
        
    if not os.path.isdir(snp_output_path):
        os.mkdir(snp_output_path)    
        
    core_gene_name_list = extract_core_gene(df_genome_info, args.output_dir, core_gene_file_path, args.quiet)
    
    merge_gap_for_msa(df_genome_info, [gn for gn in df_genome_info.index.tolist() if df_genome_info.loc[gn, 'type'] == 'complete'][0], core_gene_name_list, core_gene_file_path, delta_output_path, snp_output_path, args.quiet)    
    tree_all = generate_ml_tree('out', core_gene_file_path, args.quiet)        
    tree_final = get_final_tree(tree_all, df_genome_info, core_gene_file_path, "reuteri_core_final.nwk")

    subgroup_dict, uid_type_dict, df_snp = cluster_strain_types(df_genome_info, tree_final, args.snv_threshold, core_gene_file_path, delta_output_path, snp_output_path, args.quiet)
    
    annotate_internal_tree("reuteri_core_final.nwk", subgroup_dict, core_gene_file_path)
    tree_rep_internal = Tree("{0}/reuteri_core_rep_internal.nwk".format(core_gene_file_path), format = 1)
    
    kraken_fasta_file_path = "{0}/fasta".format(args.db)
    kraken_taxonomy_file_path = "{0}/taxonomy".format(args.db)
    
    if not os.path.isdir(args.db):
        os.mkdir(args.db)
        
    if not os.path.isdir(kraken_fasta_file_path):
        os.mkdir(kraken_fasta_file_path)
        
    if not os.path.isdir(kraken_taxonomy_file_path):
        os.mkdir(kraken_taxonomy_file_path)
    
    name_id_dict = generate_kraken_taxonomy_file(tree_rep_internal, kraken_taxonomy_file_path)    
    generate_kraken_fasta_file(subgroup_dict, name_id_dict, core_gene_file_path, args.db, kraken_fasta_file_path)    
    build_kraken_bracken_db(subgroup_dict, args.db, kraken_fasta_file_path, args.quiet)   
