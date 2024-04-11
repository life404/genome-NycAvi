#! /usr/bin/python3

from pathlib import Path
from utils import cmd_run_multiple
#from utils import make_parse
from Bio import SeqIO
import argparse
import shutil


def prank_MSA(input, output, threads):
    fasta_lst = list(Path(input).glob("*.faa"))

    root_path = Path(output)
    root_path.mkdir(exist_ok=True, parents=True)
    output = root_path / "prank"
    output.mkdir(exist_ok=True, parents=True)

    prank_cmds = [
        f'~/TOOLS/prank/bin/prank +F -d={fasta} -o={output}/{fasta.name.replace(".faa", ".prank")}'
        for fasta in fasta_lst
    ]

    check = cmd_run_multiple(prank_cmds, root_path, check="prank", num=threads)


def pal2nal(output, threads, check, mrna):
    root_path = Path(output)

    output = root_path / "pal2nal"
    output.mkdir(exist_ok=True, parents=True)

    prank_results = list((root_path / "prank").glob("*.best.fas"))
    mrna_lst = [Path(mrna) / f'{i.name.replace(".prank.best.fas", "")}.fna' for i in prank_results]

    pal2nal_cmds = []
    if len(prank_results) == len(mrna_lst):
        for i in range(0, len(prank_results)):
            pal2nal_cmds.append(
                f"~/TOOLS/pal2nal.v14/pal2nal.pl {prank_results[i]} {mrna_lst[i]} -output fasta > {output}/{prank_results[i].name.replace('.prank.best.fas', '')}.pal2nal"
            )
    else:
        print(
            f"There are something wrong, some files of results of trimAl is not in mRNA input directory"
        )
    check = cmd_run_multiple(pal2nal_cmds, root_path, check="pal2nal", num=threads)
    return check


def iqtree_gene(output, threads):
    root_path = Path(output)

    output = root_path / "iqtree_gene"
    output.mkdir(exist_ok=True, parents=True)

    # modified the fasta id in pal2nal sequences files
    iqtree_input = output / "input"
    iqtree_input.mkdir(exist_ok=True, parents=True)
    pal2nal_results = list((root_path / "pal2nal").glob("*.pal2nal"))
    for pal2nal_fa in pal2nal_results:
        content = []
        for record in SeqIO.parse(pal2nal_fa, "fasta"):
            record.id = record.id.split("|")[1]
            content.append(record)
        SeqIO.write(content, f"{iqtree_input}/{pal2nal_fa.stem}.fa", "fasta")

    iqtree_out = output / "output"
    iqtree_out.mkdir(exist_ok=True, parents=True)

    iqtree_cmds = [
        f"~/TOOLS/iqtree-2.1.3-Linux/bin/iqtree2 -s {i} -B 1000 --boot-trees --prefix {iqtree_out}/{i.stem} --quiet"
        for i in list(iqtree_input.glob("*.fa"))
    ]
    check = cmd_run_multiple(iqtree_cmds, root_path, check="iqtree", num=threads)

    return check


def ASTRAL(output):
    root_path = Path(output)

    output = root_path / "ASTRAL"
    output.mkdir(exist_ok=True, parents=True)

    astral_input = output / "input"
    astral_input.mkdir(exist_ok=True, parents=True)
    ml_best = open(f"{astral_input}/ml_best", "w")
    for contree in list((root_path / "iqtree_gene" / "output").glob("*.contree")):
        with open(contree, "r") as file:
            for line in file:
                ml_best.write(line)
    ml_best.close()

    ml_boot = open(f"{astral_input}/ml_boot", "w")
    ml_boot.writelines(
        [
            f"{i}\n"
            for i in list((root_path / "iqtree_gene" / "output").glob("*.ufboot"))
        ]
    )
    ml_boot.close()

    # cmd1: astral with ufboot files; cmd2: astral without ufboot files
    cmd1_output = output / "BEST_BOOTSTRAP"
    cmd1_output.mkdir(exist_ok=True, parents=True)
    cmd2_output = output / "BEST_ONLY"
    cmd2_output.mkdir(exist_ok=True, parents=True)

    astral_cmd1 = f"""
    java -D"java.library.path=/home/panda2bat/TOOLS/ASTRAL/lib/" -jar /home/panda2bat/TOOLS/ASTRAL/astral.5.15.5.jar -C -i {astral_input}/ml_best -b {astral_input}/ml_boot -o {cmd1_output}/astral_boot.tree > {cmd1_output}/astral_boot.log
    """
    astral_cmd2 = f"""
    java -D"java.library.path=/home/panda2bat/TOOLS/ASTRAL/lib/" -jar /home/panda2bat/TOOLS/ASTRAL/astral.5.15.5.jar -C -i {astral_input}/ml_best -b {astral_input}/ml_boot -o {cmd2_output}/astral_pp.tree > {cmd2_output}/astral_pp.log
    """
    check = cmd_run_multiple([astral_cmd1, astral_cmd2], root_path, check="astral")


def trimAl(path, ou_path, threads):
    input_path = Path(path)
    output_path = Path(ou_path)
    output_path.mkdir(exist_ok=True, parents=True)

    input_files = list(input_path.glob("*.fas"))
    trimAl_cmds = [
        f"""/home/panda2bat/TOOLS/trimAl/source/trimal \
            -in {input_file} \
            -out {output_path}/{f"{input_file.name.replace('.prank.best.fas', '')}.trimal"} \
            -nogaps"""
        for input_file in input_files
    ]

    check = cmd_run_multiple(trimAl_cmds, ou_path = output_path.parent, num=threads, check='trimAl')

    return check

def concatenated(path):
    input_path = Path(path)
    input_files = list(input_path.glob('*.trimal'))
    seq_dir = {}
    for input_file in input_files:
        records = SeqIO.parse(input_file, 'fasta')
        for record in records:
            species = record.id.split('|')[1]
            if species in seq_dir.keys():
                seq_dir[species].append(str(record.seq))
            else:
                seq_dir[species] = []
                seq_dir[species].append(str(record.seq))
    
    output_file = input_path / 'concatenated.fa'
    output_open = open(output_file, 'w')
    for key, value in seq_dir.items():
        output_open.write(f">{key}\n{''.join(value)}\n")
    output_open.close()
    return output_file
    
def estimate_length(output):
    species_tree = Path('/home/panda2bat/Avivorous_bat/output/12_evolution-species_tree/ASTRAL/BEST_BOOTSTRAP/species.tree')
    root_path = Path(output)
    trimal_path = root_path / "trimAl"
    trimal_path.mkdir(exist_ok=True, parents=True)

    align_output = root_path / "prank"
    check = trimAl(align_output, trimal_path, 128)
    
    output_path = root_path / 'ASTRAL' / 'estimated_length'
    output_path.mkdir(exist_ok=True, parents=True)
    
    if check.exists():
        concatenated_fa = concatenated(trimal_path)
        iqtree_cmd = f"""
            /home/panda2bat/TOOLS/iqtree-1.6.12-Linux/bin/iqtree -s {concatenated_fa} \
                -te {species_tree} -keep-ident -pre {output_path}/'iqtree_estimated' -nt 240 -safe
        """
        iqtree2_cmd = f"""
            ~/TOOLS/iqtree-2.1.3-Linux/bin/iqtree2 -s {concatenated_fa} \
                -t {species_tree} -m TEST --mset raxml --tree-fix --safe -T 240 -pre {output_path}/'iqtree_estimated' 
        """
        check = cmd_run_multiple(iqtree2_cmd, ou_path=output_path, check='iqtree_estimated', verbose=True)

def concatenated2(path, output):
    input_path = Path(path)
    output_path = Path(output)
    input_files = list(input_path.glob('*.pal2nal'))
    seq_dir = {}
    for input_file in input_files:
        records = SeqIO.parse(input_file, 'fasta')
        for record in records:
            species = record.id.split('|')[1]
            if species in seq_dir.keys():
                seq_dir[species].append(str(record.seq))
            else:
                seq_dir[species] = []
                seq_dir[species].append(str(record.seq))
    
    output_file = output_path / 'concatenated.fa'
    output_open = open(output_file, 'w')
    for key, value in seq_dir.items():
        output_open.write(f">{key}\n{''.join(value)[2::3]}\n")
    output_open.close()
    return output_file
    
#def divergence_time(output, tree):
#    root_path = Path(output)
#    
#    output = root_path / 'mcmctree'
#    output.mkdir(exist_ok=True, parents=True)
#    
#    mcmc_input = output / 'input'
#    mcmc_input.mkdir(exist_ok=True, parents=True)
#    
#    # generate 3rd sites file
#    pal2nal_path = root_path / 'pal2nal'
#    concatenated_fa = concatenated2(pal2nal_path, mcmc_input)
#    # copy tree file into input directory of MCMC
#    
#    # run baseml first
#    baseml_tmp = open(Path(__file__).absolute().parent / '.mcmctree_1_baseml.ctl', 'r').read()
#    baseml_tmp = baseml_tmp.replace('$seqfile$', str(concatenated_fa.absolute()))
#    baseml_tmp = baseml_tmp.replace('$treefile$', str(Path(tree).absolute()))
#    baseml_path = output / 'baseml'
#    baseml_path.mkdir(exist_ok=True, parents=True)
#    config_path = open((baseml_path / 'baseml.ctl'), 'w')
#    config_path.writelines(baseml_tmp)
#    config_path.close()
#    
#    baseml_cmd = f'/home/panda2bat/TOOLS/paml4.9j/bin/baseml'
#    check = cmd_run_multiple(baseml_cmd, ou_path=baseml_path, check='baseml', verbose=True)
#    
    
def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', '--input', dest = 'input', type = str, help = 'The input directory')
    parse.add_argument('-o', '--output', dest = 'output', type = str, help = 'The output directory')
    parse.add_argument('--threads', dest = 'threads', type = int, help = 'The number of threads')
    parse.add_argument('--step', dest = 'step', type = int, help = 'Each number corresponds to a specific step, 1) prank MSA; 2) pal2nal; 3) iqtree with every gene; 4) ASTRALIII; 5) estimated branch length with concatenated genes')
    parse.add_argument('-m', dest = 'mrna', type = str, help = 'The directory contains corresponding mRNA fasta files')
    parse.add_argument('-t', dest='tree', type=str, help='The labeled tree file used in mcmctree')
    args = parse.parse_args()
    return(args)

def main():
    args = make_parse()
    step = args.step

    if step == 1:
        prank_MSA(args.input, args.output, threads = 240)
    elif step == 2:
        mrna = args.mrna
        check = pal2nal(args.output, threads=240, check="pal2nal", mrna=mrna)
    elif step == 3:
        iqtree_gene(args.output, threads=200)
    elif step == 4:
        ASTRAL(args.output)
    elif step==5:
        estimate_length(args.output)
#    elif step==6:
#        divergence_time(args.output, args.tree)
    else:
        print("PLEASE INPUT THE CORRECT STEP NUM")


if __name__ == "__main__":
    main()
