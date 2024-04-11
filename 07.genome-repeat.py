#!/usr/bin/env python3

import argparse
from pathlib import Path
import subprocess
import time 
from utils import *

def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument('--genome', dest = 'genome', help = 'The finally genome assembly')
    parse.add_argument('-o', dest = 'ou_path', help = 'The path of output directory')
    parse.add_argument('--cds', dest = 'cds', help = 'The cds fasta of its close relative specie')
    parse.add_argument('--curatedlib', dest = 'curated', help = 'The curated lib file')
    parse.add_argument('--step', dest = 'step', type = int, help = 'Each number corresponds to a specific step, 1:EDTA; 2:softmask;')
    args = parse.parse_args()
    return(args)

#def cmd_run(command, ou_path):
#    run = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, bufsize = 1, cwd = ou_path, text = True)
#    
#    while run.poll() is None:
#        line = run.stdout.readline()
#        run.stdout.flush()
#        print(line)
#    return(run.wait())

def EDTA(genome, ou_path, cds, curatedlib):
    ou_path = Path(ou_path)
    obtain_cds_commands = f'~/TOOLS/seqkit grep -n -r -p "mRNA$" {cds} -o {ou_path}/cds.fa'
    check = cmd_run_multiple(obtain_cds_commands, ou_path, check = 'prepare_cds')
    
    if check.exists():
        EDTA_command = f'conda run -n EDTA EDTA.pl --genome {genome} --species others --step all --overwrite 0 --cds {ou_path}/cds.fa --curatedlib {curatedlib} --anno 1 --threads 64'
        check = cmd_run_multiple(EDTA_command, ou_path, check = 'EDTA_run')
   
def EDTA_soft(genome, ou_path):
    rmout_path = list(Path(ou_path).rglob('*.RM.out'))[0]
    softmask_path = Path(ou_path) / f'{Path(genome).name}.mod.SOFT.masked'
    # The softmask genome are produced using parameter '-maxdiv 30 -minscore 1000 -minlen 1000 -hardmask 0 misschar N'. The masked regions in sofemasked genome should be same as the hard masked regions in *.MAKER.masked. The `-minlen 1000` parameter contains the short repeat sequences, which may localed in genetic regions.
    ou_path = Path(ou_path)
    
    soft_command = f'perl /home/panda2bat/TOOLS/miniconda3/envs/EDTA/share/EDTA/util/make_masked.pl -genome {genome} -rmout {rmout_path} -maxdiv 30 -minscore 1000 -minlen 1000 -hardmask 0 -misschar N -t 128'
    cmd_run_multiple(soft_command, ou_path)
    
    Path(f'{genome}.new.masked').replace(softmask_path)
    

            
    
def main():
    args = make_parse()
    genome = args.genome
    ou_path = args.ou_path
    step = args.step
    
    Path(ou_path).mkdir(exist_ok=True, parents=True)
    if step == 1:
        cds = args.cds
        curatedlib = args.curated
        EDTA(genome, ou_path, cds, curatedlib)
    if step == 2:
        EDTA_soft(genome, ou_path)
    else:
        print('Please input the correct step number')
        exit(1)
        

if __name__ == '__main__':
    main()
    

