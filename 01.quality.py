#! /usr/bin/env python3
# Sequence filter
# The sequence with low quality base more than 20% or 'N' base more than 10% will be discarded
# A base with Pread <= 20 will be treaded as a low quality base

import argparse
from pathlib import Path


def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument(
        "--flist",
        type=str,
        help="The file contains names of all input files. The file contains two columns, the delimter is tab. The first column is the name of first file of PE paried, the second is the second of PE paired. \n The file should not contain the paths of any input files, the path of input file can be set uing parameter `--in`",
        required=True,
        dest="flist"
    )
    parse.add_argument(
        "--in",
        type=str,
        help="The path of directory, which contains all input files",
        dest="input"
    )
    parse.add_argument(
        "--ou",
        type=str,
        help="The ouput path, if don't set, the output path is input path",
        dest="ouput"
    )
    parse.add_argument("--na", type=str, help="The prefix of output files", dest="name")
    args = parse.parse_args()
    return args


def generate_command(f_path, o_path, o_name, f_file):
    o_num = 1
    with open (f_file) as file:
        for line in file:
            line = line.strip()
            array = line.split("\t")
            file1 = array[0]
            file2 = array[1]
            command = "fastp -i {read1} -I {read2} -o {out1} -O {out2} -q 20 -u 20 -n 15 -w 16".format(
                read1 = f_path / Path(file1),
                read2 = f_path / Path(file2),
                out1 = o_path / Path(f"{o_name}-{o_num}_1.fq.gz"),
                out2 = o_path / Path(f"{o_name}-{o_num}_2.fq.gz")
            )
            o_num += 1
            print(command)

def main():
    args = make_parse()
    f_path = Path(args.input) 
    o_path = Path(args.ouput)
    o_name = args.name 
    f_file = args.flist
    
    
    o_path.mkdir(exist_ok=True, parents=True)
    generate_command(f_path, o_path, o_name, f_file)
    
if __name__ == "__main__":
    main()