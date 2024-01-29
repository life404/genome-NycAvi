#! /usr/bin/env python3

import argparse
from pathlib import Path
from utils import *



def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument('--genome', dest ='genome', type = str, help = 'The finally assembly genome fasta file')
    parse.add_argument('--ou', dest = 'ou_path', type = str, help = 'The path of output directory')
    parse.add_argument('--protein', dest = 'protein', type = str, help = 'The protein fasta of annotation, only be used in annotaion completes analysis')
    parse.add_argument('--step', dest = 'step', type = int, help = 'The step of pargrams: 1:busco with metaeuk; 2: busco with augustus; 3: annotation completes using busco')
    args = parse.parse_args()
    return(args)

def busco_metaeuk(genome, ou_path):
    ou_path = Path(ou_path) / 'busco'
    ou_path.mkdir(exist_ok=True, parents=True)
    
    busco_command = f"/home/panda2bat/TOOLS/miniconda3/envs/busco/bin/busco -i {genome} -o metaeuk -m geno -l laurasiatheria_odb10 -c 128 --out_path {ou_path}"
    print(busco_command)
    
def busco_augustus(genome, ou_path):
    ou_path = Path(ou_path) / 'busco'
    ou_path.mkdir(exist_ok=True, parents=True)
    
    busco_command = f"/home/panda2bat/TOOLS/miniconda3/envs/busco/bin/busco -i {genome} -o augustus -m geno -l laurasiatheria_odb10 --augustus -c 128 --out_path {ou_path}"
    print(busco_command)
 
def annotation_completes(protein, ou_path):
      root = Path(ou_path)
      ou_path = root / 'busco'
      ou_path.mkdir(exist_ok=True, parents=True)
      
      #busco_command = f"conda run -n busco busco -i {protein} -o annotation_completes -m prot -l laurasiatheria_odb10 -c 128 --out_path {ou_path}"
      busco_command = f"conda run -n busco3 run_BUSCO.py -i {protein} -o annotation_completes_V3 -m prot -l /home/panda2bat/TOOLS/busco_downloads/mammalia_odb9 -c 128"
      check = cmd_run_multiple(busco_command, ou_path, check = 'busco_annotation_completes', verbose = True)
      
      
def main():
    args = make_parse()
    step = args.step
    ou_path = args.ou_path   
    
    if step == 1:
        genome = args.genome
        busco_metaeuk(genome, ou_path)
    elif step == 2:
        genome = args.genome
        busco_augustus(genome, ou_path)
    elif step == 3:
        protein = args.protein
        annotation_completes(protein, ou_path)
    else:
        print('Please choose correct step')
        exit(1)
        
if __name__ == "__main__":
    main()
    
    

    
    
    