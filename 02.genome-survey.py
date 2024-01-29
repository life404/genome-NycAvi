#! /usr/bin/env python3

import argparse
import subprocess
from pathlib import Path

def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument("--flist", dest = "flist", type = str, help = "The file contains all input, each row was divided two columns by 'tab', the first column was the first file of a paired-end fastq file, and the second row is the second")
    parse.add_argument("--in", dest = "input", type = str, help = "the path of the directory, which contains all fastq files")
    parse.add_argument("--ou", dest = "output", type = str, help = "the path of output directory of resutls")
    parse.add_argument("--na", dest = "name", type = str, help = "the prefix of results files")
    parse.add_argument("--range", dest = "krange", type = str, help = "the range of kmer sizes, the format of kmer range is start:end:seq, the programe will walk all kmer values.\nThe default value of kmer range is 21 to 21, the seq is 1", default = "21:21:1")
    parse.add_argument("--jpath", dest = "jpath", type = str, help = "The absoult path of jellyfish programe", default = "/home/panda2bat/TOOLS/miniconda3/bin/jellyfish")
    args = parse.parse_args()
    return args

def generage_command(in_path, ou_path, name, k_range, jellyfish_path):
    
    flist = " ".join([str(i) for i in list(Path(in_path).glob("*.fq.gz"))])
    
    for kmer in range(int(k_range.split(":")[0]), int(k_range.split(":")[1]) + 1, int(k_range.split(":")[2])):
        results = Path(ou_path) / Path(f"{name}.{kmer}kmer.jf")
        results_2 = Path(ou_path) / Path(f"{name}.{kmer}kmer.histo")
        jellyfish_count = f"{jellyfish_path} count -C -m {kmer} -s 100G -t 128 -o {results} <(zcat {flist})"
        jellyfish_histo = f"{jellyfish_path} histo -t 64 {results} > {results_2}"
        print(f"{jellyfish_count}\n{jellyfish_histo}\n")
        

def main():
    args = make_parse()
    in_path = args.input 
    ou_path = args.output
    name = args.name 
    k_range = args.krange
    jellyfish_path = args.jpath
    
    Path(ou_path).mkdir(exist_ok = True, parents= True)
    
    generage_command(in_path, ou_path, name, k_range, jellyfish_path)
    
if __name__ == "__main__":
    main()
    
    
        
        
    
    
    