#! /usr/bin/python3 

"""
这个脚本用于从NCBI下载的蛋白质序列提取最长的蛋白质序列和转录本序列用于orthofinder以及系统发育的构建
NCBI下载请使使用datasets

datasets download genome accession GCF_000001405.40 --include gff3,cds,protein

"""

import gffutils
from pathlib import Path
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd 
import argparse
import subprocess


def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument('-g', '--gff', dest='gff', help='The NCBI gff3 file')
    #parse.add_argument('-c', '--cds', dest='cds', help='The raw cds file')
    #parse.add_argument('-p', '--pep', dest='pep', help='The raw pep file')
    parse.add_argument('-f', '--fna', dest='fna', help='The geneome file')
    parse.add_argument('--prefix', dest='prefix', default='the_longest', help='The prefix of the longest files, such as: prefix.gff3, prefix.cds, prefix.pep')
    args = parse.parse_args()
    return args


def parse_gff3(args):
    
    if not Path('gff.db').exists():
        db = gffutils.create_db(str(Path(args.gff).resolve()), dbfn = 'gff.db', force=True, keep_order=True, id_spec={'gene':'ID', 'mRNA':'ID'})
    else:
        db = gffutils.FeatureDB('gff.db')
    return db

#def parse_fasta(args):
#    cds_lst = SeqIO.parse(Path(args.cds).resolve(), 'fasta')
#    pep_lst = SeqIO.parse(Path(args.pep).resolve(), 'fasta')
#
#    # 只保留蛋白质编码基因
#    cds_records = {re.search(r'protein_id=[\w\d.]*', record.description).group().split('=')[1]:record for record in cds_lst if 'protein_id' in record.description}
#    #for record in cds_lst:
#    #    print(record.description)
#    #    re.search(r'protein_id=[\w\d.]*', record.description).group()
#        
#    pep_records = {record.id:record for record in pep_lst}
#
#    return cds_records, pep_records

def keep_longest(db, args):

    out_gff = Path(f"{args.prefix}.gff3").resolve()
    out_cds = Path(f"{args.prefix}.cds").resolve()
    out_pep = Path(f"{args.prefix}.pep").resolve()
    # 保留蛋白质编码基因
    genes = [gene for gene in db.features_of_type('gene') if gene.attributes['gene_biotype'][0]=='protein_coding']

    longest_cds = []
    longest_pep = []
    gff_buff = open(out_gff, 'w')

    num = 1
    for gene in genes:

        mRNA = [i for i in db.children(gene.id, featuretype = 'mRNA')]
        max_mRNA = max(mRNA, key = lambda x: db.children_bp(x.id, child_featuretype='CDS'))

        first_CDS = next(db.children(max_mRNA.id, featuretype='CDS'))
        protein_id = first_CDS.attributes['protein_id'][0]
        
        # 获取基因名称，并对输出的CDS和蛋白质序列的名称进行修改
        gene_name = first_CDS.attributes['gene'][0]

        gff_buff.write(f"{gene}\n")
        gff_buff.write(f"{max_mRNA}\n")
        for i in db.children(max_mRNA.id, featuretype = 'exon'):
            gff_buff.write(f'{i}\n')
        for i in db.children(max_mRNA.id, featuretype = 'CDS'):
            gff_buff.write(f'{i}\n')
        print(f'Parsing {num} gene', end='\r')
        num += 1

    gff_buff.close()

    return out_gff

    # NCBI下载的蛋白质和CDS序列存在不匹配的情况，因此使用EVM脚本提取对应的pep和CDS
    #SeqIO.write(longest_cds, out_cds, 'fasta')
    #SeqIO.write(longest_pep, out_pep, 'fasta')

def obtain_pep_cds(out_gff, args):
    
    out_pep = Path(f"{args.prefix}.pep").resolve()
    out_cds = Path(f"{args.prefix}.cds").resolve()

    command_cds = f'~/TOOLS/EVidenceModeler-v2.1.0/EvmUtils/gff3_file_to_proteins.pl {out_gff} {args.fna} prot > {out_pep}'
    command_pep = f'~/TOOLS/EVidenceModeler-v2.1.0/EvmUtils/gff3_file_to_proteins.pl {out_gff} {args.fna} CDS  > {out_cds}'

    process_cds = subprocess.Popen(command_pep, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
    process_pep = subprocess.Popen(command_cds, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)

    if process_cds.wait() != 0:
        print(f"ERROR: {process_cds.communicate()[1]}")
    if process_pep.wait() != 0:
        print(f"ERROR: {process_pep.communicate()[1]}")

    return out_pep, out_cds

def remove_remainder(record, remainder):
    record_list = list()
    for start in range(0, remainder+1):
        end = None if (remainder - start) == 0 else -(remainder - start)
        record_list.append(
            record[start:end]
        )
    return record_list


def fix_cds(out_pep, out_cds):
    cds_records = {i.id:i for i in SeqIO.parse(out_cds, 'fasta')}
    pep_records = {i.id:i for i in SeqIO.parse(out_pep, 'fasta')}

    fixed = list()
    error = list()
    num = 1
    for key, value in cds_records.items():
        target = pep_records.get(key)
        remainder = len(value.seq) - len(target.seq)*3
        if remainder != 0:
            flag = False
            for i in remove_remainder(value, remainder):
                if i.translate().seq == target.seq:
                    fixed.append(
                        SeqRecord(
                            i.seq,
                            id=i.id,
                            name=i.name,
                            description=i.description
                        )
                    )
                    flag = True

            if not flag:
                error.append(value)
            
        else:
            fixed.append(value)

        print(f"Fixed {num} cds seqs", end='\r')
        num += 1

    SeqIO.write(fixed, out_cds, 'fasta')
    
    if error:
        SeqIO.write(error, 'error.cds', 'fasta')

def main():
    args = make_parse()
    db = parse_gff3(args)
    out_gff = keep_longest(db, args)
    print(f"Finsih extrac the longest isoforms")
    out_pep, out_cds = obtain_pep_cds(out_gff, args)
    print(f"Finsih obtain protein and cds sequences")

    fix_cds(out_pep, out_cds)

if __name__ == "__main__":
    main()

