#! /usr/bin/env python3

import argparse
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd

def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument(
        "--in",
        dest="input",
        type=str,
        help="the path of input directory, which contains all input files",
    )
    parse.add_argument(
        "--ou", dest="output", type=str, help="the path of output directory"
    )
    parse.add_argument(
        "--gsize",
        dest="gsize",
        type=str,
        help="the prediceted size of genome, a string ended with g, for example, 2g",
    )
    parse.add_argument(
        "--step",
        dest="step",
        type=int,
        help="the different step of programe, 1 means nextDenovo; 2 means kbmap of wtdbg2; 3 means wtdbg of wtdbg2"
    )
    args = parse.parse_args()
    return args


def nextDenovo_command(in_path, ou_path, gsize):
    ou_path = Path(ou_path) / "nextDenovo"
    ou_path.mkdir(exist_ok=True, parents=True)

    flist = [str(i) + "\n" for i in list(in_path.rglob("*.fasta.gz"))]
    fofo_path = Path(ou_path) / "fofo"
    fofo = open(fofo_path, "w")
    fofo.writelines(flist)
    fofo.close()

    next_config_tmp = open(Path(__file__).parent / ".nextDenovo.config").read()
    next_config_tmp = next_config_tmp.replace("$FOFO$", str(fofo_path))
    next_config_tmp = next_config_tmp.replace("$OUTPUT$", str(ou_path))
    next_config_tmp = next_config_tmp.replace("$GSIZE$", gsize)

    config_path = Path(ou_path) / "correct.cfg"
    config = open(config_path, "w")
    config.writelines(next_config_tmp)
    config.close()

    nextDenovo_run = f"/home/panda2bat/TOOLS/NextDenovo/nextDenovo {config_path}"
    print(nextDenovo_run)

def wtdbg_command(ou_path):
    cns_path = list(Path(ou_path).joinpath("nextDenovo").rglob("cns.fasta"))

    # make softlink of cns fasta from nextDenovo results
    ou_path = Path(ou_path) / "wtdbg"
    ou_path.mkdir(exist_ok=True, parents=True)

    for cns in cns_path:
        link_path = ou_path / (str(cns.parts[-2]) + ".fasta")
        link_path.symlink_to(cns,)
        
        
    kmp_input_list = list(ou_path.rglob("seed_cns*.fasta"))
    kmp_run_lst = []
    
    for i in kmp_input_list:
        for j in kmp_input_list:
            kmp_out = (
                Path(ou_path) / ("kbmap."
                + str(i.stem)
                + "."
                + str(j.stem)
            ))
            kmp_run_lst.append(
                f"/home/panda2bat/TOOLS/wtdbg-2.5_x64_linux/kbm2 -t 10 -c 10 -i {str(i)} -d {str(j)} -fo {str(kmp_out)}\n"
            )
    kmp_command = open(
        Path(__file__).parent / "03_genome-assembly-wtdbg2-kmap.jobs", "w"
    )
    kmp_command.writelines(kmp_run_lst)
    kmp_command.close()

    return kmp_input_list


def wtdbg_command2(ou_path, gsize):
    ou_path = Path(ou_path) / "wtdbg"
    kmp_lst = [str(i) for i in list(ou_path.glob("kbmap*"))]
    kmp_fasta_list = " -i ".join([str(i) for i in list(ou_path.glob("seed_cns*.fasta"))])

    merge_kmp_path = Path(ou_path) / "multip-node.kbmap"
    if not merge_kmp_path.exists():
        merge_kmp = subprocess.Popen(
            f'cat {" ".join(kmp_lst)} > {str(ou_path / "multip-node.kbmap")}', shell=True
        )
        print(merge_kmp.wait())
        if merge_kmp.poll() != 0:
            exit(1)

    wtdbg_run_lst = []
    for nodelen in [1536, 2048, 2304, 2560, 1024]:
        for s in [0.5]:
            for alndovetail in [-1, 4608, 9216, 256]:
                for L in [5000, 0]:
                    dir_path = (
                    Path(ou_path)
                        / f"ND0.25.NL{nodelen}.NM400.S{s}.E3.ALN{alndovetail}.L{L}"
                    )
                    dir_path.mkdir(exist_ok=True, parents=True)
                    
                    out_path = (
                        Path(ou_path)
                        / f"ND0.25.NL{nodelen}.NM400.S{s}.E3.ALN{alndovetail}.L{L}/wtdbg"
                    )
                    log_path = (
                    Path(ou_path)
                        / f"ND0.25.NL{nodelen}.NM400.S{s}.E3.ALN{alndovetail}.L{L}/wtdbg.log" 
                    )
                    wtdbg_run_lst.append(
                        f'/home/panda2bat/TOOLS/wtdbg-2.5_x64_linux/wtdbg2 -i {kmp_fasta_list} --load-alignments {merge_kmp_path} -g {gsize} -t 128 -K 0.05 -A --node-drop 0.25 --node-len {nodelen} --node-max 400 -s {s} -e 3 -L {L} --rescue-low-cov-edges --no-read-length-sort --aln-dovetail {alndovetail} -fo {out_path} &> {log_path} \n'
                    )

    wtdbg_command = open(
        Path(__file__).parent / "03_genome-assembly-wtdbg2-wtdbg.jobs", "w"
    )
    wtdbg_command.writelines(wtdbg_run_lst)
    wtdbg_command.close()
    
def obtain_longest(ou_path):
    ou_path = Path(ou_path) / 'wtdbg'
    result_dir = [dir for dir in ou_path.iterdir() if dir.is_dir()]
    
    infor = []
    for dir in result_dir:
        log_path = dir / 'wtdbg.log'
        with open (log_path) as file:
            for line in file:
                if 'Estimated' in line:
                    infor.append([int(i.split(" ")[-1]) for i in line.split(",") ])
        
    infor_data = pd.DataFrame(infor, columns = ["TOT", "CNT", "AVG", "MAX", "N50", "L50", "N90", "L90", "Min"])
    infor_data['PARAMETER'] = [str(i.name) for i in result_dir]
    
    infor_data = infor_data.sort_values(by = ['TOT', 'N50', 'CNT'], ascending = (False, False, True)).reset_index(drop=True)
    print(infor_data)
    
    print(f'the longest result of wtdbg2 using parameter {infor_data.at[0, "PARAMETER"]}')
    return(infor_data.at[0, 'PARAMETER'])

def wtdbg_cns_command(ou_path, longest_one):
    ou_path = Path(ou_path) / 'wtdbg'
    longest = ou_path / longest_one / 'wtdbg.ctg.lay.gz'
    ou_fa = ou_path / f'{longest_one}.fa'
    
    cns_command= f'/home/panda2bat/TOOLS/wtdbg-2.5_x64_linux/wtpoa-cns -t 256 -i {longest} -fo {ou_fa}'
    
    print(cns_command)
    
def main():
    args = make_parse()
    in_path = Path(args.input)
    ou_path = Path(args.output)
    gsize = args.gsize
    step = args.step

    Path(ou_path).mkdir(exist_ok=True, parents=True)
    
    if step == 1:
        nextDenovo_command(in_path, ou_path, gsize)
    elif step == 2:
        wtdbg_command(ou_path)
    elif step == 3:
        wtdbg_command2(ou_path, gsize)
    elif step == 4:
        longest_one = obtain_longest(ou_path) 
        wtdbg_cns_command(ou_path, longest_one)
    else:
        print("please input correct step number")
        exit(1)
    


if __name__ == "__main__":
    main()


