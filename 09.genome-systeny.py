#! /usr/bin/python3

import argparse
from pathlib import Path
import subprocess
from utils import *
from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import globalenv
import re

def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument(
        "--target",
        "-t",
        dest="target",
        type=str,
        help="The softmasked target genome used to perform lastz alignment. The genome files should be 2bit format",
    )
    parse.add_argument("--query", "-q", dest="query", type=str, help="The softmasked query genome, 2bit format")
    parse.add_argument("--output", "-o", dest = "ou_path", type = str, help = "The path of output directory")
    parse.add_argument("--tprefix", dest = "tprefix", required=False, type = str, help = "The name of target")
    parse.add_argument("--qprefix", dest = "qprefix", required=False, type = str, help = "The name of query")
    parse.add_argument("--step", dest = "step", type = str, help = "The steps of program, the step 1:lastz alignment;")
    parse.add_argument("--verbose", dest = "verbose", type = bool, default = True, choices=[True, False], help = "Print all (True, default) information of programe or not (False)")
    parse.add_argument("--threads", dest = "threads", type = int, default = 5, help = "The nubmer of threads")
    args = parse.parse_args()
    return(args)

def lastz_alignment(target, query, tprefix, qprefix, ou_path, verbose, threads):
    root_path = Path(ou_path) / 'lastz'
    root_path.mkdir(exist_ok=True, parents=True)
    target = Path(target)
    query = Path(query)
    all_chain = root_path / 'all.chain'
    preNet_chain = root_path / 'all.pre.chain'
   
    # Beforing run this function, please load the ucsc envirment using module 
    chains_dir = root_path / 'chains'
    chains_dir.mkdir(exist_ok=True, parents=True)
    # Initialize the R console
    rconsole = robjects.r
    rconsole_initialize = f"""
    pacman::p_load(parallel, CNEr)
    assemblyTarget <- '{target}'
    assemblyQuery  <- '{query}'
    outputDir <- '{chains_dir}'
    """ 
    rconsole(rconsole_initialize)
    
    if not (root_path / 'tmp.jobs').exists(): 
        lastz_r_aligner = f"""
        lavs <- lastz(assemblyTarget = assemblyTarget, assemblyQuery = assemblyQuery, distance = 'near', outputDir = outputDir, echoCommand = TRUE)
        """
        rconsole(lastz_r_aligner)
        lastz_r_command_lst = list(rconsole.lavs)
        lastz_r_command_file = open(root_path / 'tmp.jobs', 'w')
        for i in lastz_r_command_lst:
            lastz_r_command_file.write(i + '\n')
        lastz_r_command_file.close()
    
    lastz_r_command_lst = open(root_path / 'tmp.jobs', 'r').read()
    check = cmd_run_multiple(lastz_r_command_lst, root_path, check = 'lastz_generate_aligner', num = threads, verbose = verbose)
    
    if check.exists():
        lavs = [str(i) for i in list(chains_dir.glob('*.lav'))]
        globalenv['lavs'] = robjects.StrVector(lavs)
        lastz_r_Psl = """
        psls <- lavToPsl(lavs, removeLav = FALSE, binary = 'lavToPsl')
        """
        check = cmd_run_R(rconsole, lastz_r_Psl, root_path, 'lastz_lavToPsl')
        #print(rconsole.psls)
    if check.exists():
        psls = [str(i) for i in list(chains_dir.glob('*.psl'))]
        globalenv['psls'] = robjects.StrVector(psls)
        lastz_r_chaining = f"""
        chains <- axtChain(psls, assemblyTarget  = assemblyTarget, assemblyQuery = assemblyQuery, distance = "near", removePsl = FALSE, binary = "axtChain")
        allChain <- chainMergeSort(chains, assemblyTarget = assemblyTarget, assemblyQuery = assemblyQuery, allChain = '{all_chain}', removeChains = FALSE, binary = "chainMergeSort")
        """ 
        check = cmd_run_R(rconsole, lastz_r_chaining, root_path, 'lastz_chaining')
    if check.exists():
        getSize_command = []
        getSize_command.append(f"twoBitInfo {Path(target)} .{tprefix}.info")
        getSize_command.append(f"twoBitInfo {Path(query)} .{qprefix}.info")
        check = cmd_run_multiple(getSize_command, root_path, check = 'getSize')
    if check.exists():
        chainPreNet_command = f"chainPreNet {all_chain} .{tprefix}.info .{qprefix}info {preNet_chain}"
        check = cmd_run_multiple(chainPreNet_command, root_path, check = 'lastz_ChainPreNet', verbose = verbose)
    if check.exists():
        lastz_ChainNet = f"chainNet {preNet_chain} .{tprefix}.info .{qprefix}.info stdout /dev/null | netSyntenic stdin noClass.net"
        check = cmd_run_multiple(lastz_ChainNet, root_path, check = 'lastz_chainNet', verbose = verbose)
    if check.exists():
        lastz_netToAxt = f"netToAxt noClass.net all.pre.chain {target} {query} all.axt"
        check = cmd_run_multiple(lastz_netToAxt, root_path, check = 'lastz_netToAxt', verbose = verbose)
    if check.exists():
        lastz_axtSort = f"axtSort all.axt all.sort.axt"
        check = cmd_run_multiple(lastz_axtSort, root_path, check = 'lastz_axtSort')
    if check.exists():
        lastz_axtToMaf = f"axtToMaf all.sort.axt .{tprefix}.info .{qprefix}.info {tprefix}_{qprefix}.maf -tPrefix={tprefix}. -qPrefix={qprefix}."
        check = cmd_run_multiple(lastz_axtToMaf, root_path, check = 'lastz_axtToMaf', verbose = verbose)
    
    print("╰(○'◡'○)╮: ALL LASTZ ALIGNMENT COMMANDS HAVE COMPLETED")

def lastz_alignment_update(target, query, tprefix, qprefix, ou_path, verbose, threads):
    root_path = Path(ou_path)
    ou_path = root_path / 'lastz_updated'
    ou_path.mkdir(exist_ok=True, parents=True)
    
    chain = ou_path / f'{tprefix}-{qprefix}.chainMergeSort.chain'
    Filler_chain = ou_path / f'{tprefix}-{qprefix}.RepeatFiller.chain'
    Cleaner_chain = ou_path / f'{tprefix}-{qprefix}.chainCleaner.chain'
    
    
    # using CNEr package to generate the LASTZ axt commands
    check = ou_path / 'CNEr.jobs.generate.ok'
    axtDirs = ou_path / 'axtDirs'
    axtDirs.mkdir(exist_ok=True, parents=True)
    
    if not check.exists():
        rconsole = robjects.r
        rconsole_initialize = f"""
        pacman::p_load(CNEr)
        assemblyTarget <- '{target}'
        assemblyQuery  <- '{query}'
        outputDir <- '{axtDirs}'
        """
        rconsole(rconsole_initialize)
        cner_run = f"""
        lavs <- lastz(assemblyTarget = assemblyTarget, assemblyQuery = assemblyQuery, distance = 'near', outputDir = '{axtDirs}', echoCommand = TRUE)
        """ 
        rconsole(cner_run)
        cner_command_lst = list(rconsole.lavs)
        
        cner_command_file = open(ou_path / 'CNEr.jobs', 'w')
        for i in cner_command_lst:
            i = re.sub(r'Q=.*lastzMatrix', f"Q={Path(__file__).absolute().parent.joinpath('.lastz.Matrix')}", i)
            i = i.replace('lav', 'axt')
            cner_command_file.write(i + '\n')
        cner_command_file.close()
        check.touch()
    if check.exists():
        cner_command_lst = open(ou_path / 'CNEr.jobs', 'r').readlines()
        check = cmd_run_multiple(cner_command_lst, ou_path, check = 'lastz.generate.axt', num = threads, verbose = verbose)
    if check.exists():
        axtFile_lst = [str(i) for i in axtDirs.glob('*.axt')]
        axtFile_commands = [f"axtChain -minScore=3000 -linearGap=medium {i} {target} {query} {str(i).replace('.axt', '.chain')}" for i in axtFile_lst]
        check = cmd_run_multiple(axtFile_commands, ou_path = ou_path, check = 'lastz.axtChain', num = threads, verbose = verbose)
    if check.exists():
        #chainMergeSort_commands = f"chainMergeSort {axtDirs}/*.chain > {chain}"
        # Using parameter 'inputList' to avoid the error `Argument list too long`
        chain_path_lst = [str(i)+'\n' for i in axtDirs.glob('*.chain')]
        chain_path_file = open(ou_path / 'chain_path_file', 'w')
        chain_path_file.writelines(chain_path_lst)
        chain_path_file.close()
        chainMergeSort_commands = f"chainMergeSort -inputList=chain_path_file > {chain}"
        check = cmd_run_multiple(chainMergeSort_commands, ou_path, check = 'chainMergeSort', verbose = verbose)
    if check.exists():
        RepeatFill_commands = f"RepeatFiller.py -c {chain} -T2 {target} -Q2 {query} --output {Filler_chain}"
        check = cmd_run_multiple(RepeatFill_commands, ou_path, check = 'RepeatFiller', verbose = verbose)
    if check.exists():
        TwoBitInfo_commands = [f"twoBitInfo {i} stdout | sort -k2rn > {str(i).replace('2bit', 'info')} " for i in [target, query]]
        cmd_run_multiple(TwoBitInfo_commands, ou_path, check = 'faTwoBitInfo', verbose = verbose)
        chainCleaner_commands = f"~/TOOLS/UCSCtools/chainCleaner {Filler_chain} -tSizes={str(target).replace('2bit', 'info')} -qSizes={str(query).replace('2bit', 'info')} {target} {query} {Cleaner_chain} {str(Cleaner_chain).replace('.chain', '.bed')} -linearGap=medium -doPairs  -debug"
        check = cmd_run_multiple(chainCleaner_commands, ou_path, check = 'chainCleaner', verbose = verbose)
    if check.exists():
        chainNet_commands = f"chainNet {Cleaner_chain} -minSpace=1 {str(target).replace('2bit', 'info')} {str(query).replace('2bit', 'info')} stdout /dev/null | netSyntenic stdin noClass.net"
        check = cmd_run_multiple(chainNet_commands, ou_path, check='chainNet')
    if check.exists():
        netToAxt = f"netToAxt noClass.net {Cleaner_chain} {target} {query} stdout | axtSort stdin {tprefix}-{qprefix}.axt"
        check = cmd_run_multiple(netToAxt, ou_path, check = 'netToAxt', verbose = verbose)
    if check.exists():
        axtToMaf = f"axtToMaf {tprefix}-{qprefix}.axt {str(target).replace('2bit', 'info')} {str(query).replace('2bit', 'info')} {tprefix}-{qprefix}.maf -tPrefix={tprefix}. -qPrefix={qprefix}."
        check = cmd_run_multiple(axtToMaf, ou_path, check = 'axtToMaf', verbose = verbose)
    
    print("╰(○'◡'○)╮: ALL LASTZ ALIGNMENT COMMANDS HAVE COMPLETED")   

def lastal_alignment(ou_path, target, query, tprefix, qprefix, verbose, threads):
    root_path = Path(ou_path)
    ou_path = root_path / 'lastal'
    ou_path.mkdir(exist_ok=True, parents=True)
    lastdb_command = f"conda run -n last lastdb -P 64 -uNEAR {tprefix}db {target}"
    check = cmd_run_multiple(lastdb_command, ou_path, check = 'lastdb', verbose = verbose)
    
    # peforming lastal train
    if check.exists():
        db_path = root_path / 'lastal' / f'{tprefix}db'
        ou_path = ou_path / f'{tprefix}-{qprefix}'
        ou_path.mkdir(exist_ok=True, parents=True)
        last_train = f"conda run -n last last-train -P64 --revsym -E0.05 -C2 {db_path} {query}> {tprefix}-{qprefix}.train"
        check = cmd_run_multiple(last_train, ou_path, check = 'last-train', verbose = verbose)
    if check.exists():
        lastal_command = f'parallel-fasta -j {threads} "lastal -p {tprefix}-{qprefix}.train {db_path} | last-split -fMAF+" < {query} > {tprefix}-{qprefix}.maf'
        check = cmd_run_multiple(lastal_command, ou_path, check = 'lastal', verbose = verbose)
    
    if check.exists():
       print("╰(○'◡'○)╮: ALL LAST ALIGNMENT COMMANDS HAVE COMPLETED") 
        
def main():
    args = make_parse()
    target = args.target
    query = args.query
    if not args.tprefix:
        tprefix = Path(target).stem
    else:
        tprefix = args.tprefix
    if not args.qprefix:
        qprefix = Path(query).stem
    else:
        qprefix = args.qprefix
    ou_path = args.ou_path
    step = args.step
    verbose = args.verbose
    threads = args.threads
    
    lastz_alignment_update(target, query, tprefix, qprefix, ou_path, verbose, threads)
    #lastal_alignment(ou_path, target, query, tprefix, qprefix, verbose, threads)

if __name__ == "__main__":
    main() 
    
    
    
    