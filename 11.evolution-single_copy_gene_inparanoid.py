#! /usr/bin/python3

import argparse
from utils import cmd_run_multiple
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO
import concurrent.futures as cf


def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', dest = 'input', type = str, help = 'The input directory')
    parse.add_argument('-o', dest = 'output', type = str, help = 'The output directory')
    parse.add_argument('-m', dest = 'mrna', type = str, help = 'The corresponding mRNA directory')
    parse.add_argument('--step', dest = 'step', type = int, help = 'Each number corresponds to a specific step, 1:inparanoid; 2: parse the results of inparanoid')
    args = parse.parse_args()
    return(args)

def run_inparanoid(input, output):
    input = Path(input).absolute()
    if not input.exists():
        print(f'The input directory {input} is not exist.')
        exit(1)
    output = Path(output).absolute()
    output.mkdir(exist_ok=True, parents=True)
    
    inparanoid_command = f"singularity run --bind {input}:/input/ --bind {output}:/output/ /home/panda2bat/TOOLS/inparanoid_docker/inparanoid_autorun.sif"
    check = cmd_run_multiple(inparanoid_command, ou_path = output, check = 'inparanoid', verbose = True)
    
def func(x):
    # function in apply
    ref_id = x.iloc[0, 4].split('|')[0]
    qry_id = x.iloc[1, 4].split('|')[0]
    return pd.DataFrame({x.iloc[0, 2]: [ref_id],x.iloc[1, 2]: [qry_id]}, index = None)

def parse(key, output):
    file_parse = []
    for tmp in list(output.glob(f'SQLtable.{key}*')):
        tmp_results = pd.read_csv(tmp, header = None, index_col = None, names = ['Group-id', 'Score', 'Species', 'Confidence', 'Protein'], sep = '\t')
        tmp_results = tmp_results.groupby(by = 'Group-id').filter(lambda x: len(x) == 2).groupby(by = 'Group-id').apply(func)
        file_parse.append(tmp_results)
    return key, file_parse
    
def parse_results(output):
    output = Path(output) 
    results = list(output.glob('SQLtable*'))

    # The results were be sorted
    results_species = [i.name.split('-')[0].replace('SQLtable.', '') for i in results]
    species_count = {species: results_species.count(species) for species in set(results_species)}
    species_sorted = sorted(species_count.items(), key = lambda x: x[1], reverse=True)
    # Following the descending order of number of files to trnaforme those results and merge them
    results_parse = {}
    with cf.ProcessPoolExecutor(max_workers=32) as e:
            process_lst = [e.submit(parse, key[0], output) for key in species_sorted]
            i = 1   
            for process in cf.as_completed(process_lst):
                key, file_parse = process.result()
                results_parse[key] = file_parse
                print(f'\rTotal {len(species_sorted)} jobs, have finished {i}', end = '')
                i += 1
    
    results_merge = []
    for key in species_sorted:
        key = key[0]
        if len(results_parse[key]) > 1:
            final_tmp = results_parse[key][0]
            for tmp in results_parse[key][1:]:
                final_tmp = pd.merge(final_tmp, tmp, on = [key], how = 'inner')
        else:
            final_tmp = results_parse[key][0]
        results_merge.append(final_tmp)
    
    final_results = results_merge[0]
    for tmp in results_merge[1:]:
        merge_index = list(set(tmp.columns).intersection(set(final_results)))
        final_results = pd.merge(final_results, tmp, how = 'inner', on = merge_index)
    final_results.index = [f'group{i:04d}' for i in final_results.index]
    final_results.to_csv(output / f'single_copy_gene_of_{len(species_sorted)}.csv', sep='\t')
                
    return(final_results)

def single_copy_fasta(final_results, input, output, mrna):
    root_path = Path(output)
    protein_files = Path(input).glob('*.faa')
    protein_contents = {}
    for fasta_file in protein_files:
        records = SeqIO.parse(fasta_file, 'fasta')
        protein_contents[fasta_file.name] = {}
        for record in records:
            protein_contents[fasta_file.name][record.id.split('|')[0]] = record
    
    mrna_files = Path(mrna).glob('*.fna')
    mrna_contents = {}
    for fasta_file in mrna_files:
        records = SeqIO.parse(fasta_file, 'fasta')
        mrna_contents[fasta_file.name] = {}
        for record in records:
            mrna_contents[fasta_file.name][record.id.split('|')[0]] = record
    
    output = root_path / 'single_copy_protein'
    output.mkdir(exist_ok=True, parents=True)
    protein_seq_lst = {}
    for index, data in final_results.iterrows():
        for label in data.index:
            if index in protein_seq_lst.keys():
                protein_seq_lst[index].append(protein_contents[label][data[label]])
            else:
                protein_seq_lst[index] = []
                protein_seq_lst[index].append(protein_contents[label][data[label]])
    for key, value in protein_seq_lst.items():
        SeqIO.write(value, f'{output}/{key}.faa', 'fasta')
        
    output = root_path / 'single_copy_cds'
    output.mkdir(exist_ok=True, parents=True)
    mrna_seq_lst = {}
    for index, data in final_results.iterrows():
        for label in data.index:
            if index in mrna_seq_lst.keys():
                mrna_seq_lst[index].append(mrna_contents[label.replace('.faa', '.fna')][data[label]])
            else:
                mrna_seq_lst[index] = []
                mrna_seq_lst[index].append(mrna_contents[label.replace('.faa', '.fna')][data[label]])
    for key, value in mrna_seq_lst.items():
        SeqIO.write(value, f'{output}/{key}.fna', 'fasta') 
    return(protein_seq_lst, mrna_seq_lst)
    

def main():
    args = make_parse()
    step = args.step
    input = args.input
    output = args.output
    if step == 1:
        run_inparanoid(args.input, args.output)
    elif step == 2:
        mrna = args.mrna
        final_results = parse_results(output)
        single_copy_fasta(final_results, input, output, mrna)
    else:
        print('Please set the correct step number')
        exit(1)
        
if __name__ == '__main__':
    main()
    
