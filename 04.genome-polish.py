#! /usr/bin/env python3

import argparse
from pathlib import Path

def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument("--initial", dest = "initial", help = "the initial assembly fasta from wtdbg2")
    parse.add_argument("--inl", dest = "lgs_path", help = "the input path for directory, which contains long reads")
    parse.add_argument("--ins", dest = "sgs_path", help = "the input path for directory, which contains short reads")
    parse.add_argument("--ou", dest = "ou_path", help = "the output directory of nextPolish")
    args = parse.parse_args()
    return(args)

def nextPolish_command(initial, lgs_path, sgs_path, ou_path):
    config = Path(__file__).parent / '.nextPolish.config'
    config_tmp = open(config, 'r').read()
    
    sgs_lst = [str(i) + "\n" for i in list(Path(sgs_path).glob('*.fq.gz'))]
    lgs_lst = [str(i) + "\n" for i in list(Path(lgs_path).glob('*.fasta.gz'))]
    
    sgsfofo = open(Path(ou_path) / 'sgs.fofo', 'w')
    sgsfofo.writelines(sgs_lst)
    sgsfofo.close()
    
    lgsfofo = open(Path(ou_path) / 'lgs.fofo', 'w')
    lgsfofo.writelines(lgs_lst)
    lgsfofo.close()
    
    config_tmp = config_tmp.replace('$GENOME$', initial)
    config_tmp = config_tmp.replace('$OUTPUT$', ou_path)
    config_tmp = config_tmp.replace('$SGSFOFO$', str(Path(ou_path) / 'sgs.fofo'))
    config_tmp = config_tmp.replace('$LGSFOFO$', str(Path(ou_path) / 'lgs.fofo'))
    
    polish_config = open(Path(ou_path) / 'polish.config', 'w')
    polish_config.writelines(config_tmp)
    polish_config.close()
    
    command = f"/home/panda2bat/TOOLS/NextPolish/nextPolish {str(Path(ou_path) / 'polish.config')}"
    
    print("the modified config file of nextPolish is:")
    print(config_tmp)
    print("#############################################################")
    print("Run following command to preform nextPolish:")
    print("#############################################################")
    print(command)
    
def main():
    args = make_parse()
    initial = args.initial
    lgs_path = args.lgs_path
    sgs_path = args.sgs_path
    ou_path = args.ou_path
    
    Path(ou_path).mkdir(exist_ok = True, parents = True)
    
    nextPolish_command(
		initial,
		lgs_path,
		sgs_path,
		ou_path
	)
    
if __name__ == "__main__":
    main()
    
    