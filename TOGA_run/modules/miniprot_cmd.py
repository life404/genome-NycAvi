#! /usr/bin/python3

from pathlib import Path
from modules.parallel_execute import parallel_run_cmds
from configparser import ConfigParser
import sys

config = ConfigParser()
config.read(Path(__file__).absolute().parent / "config.ini")
miniprot = Path(config["Software"]["miniprot"])


def miniprot_index(target, run_dir, check_path=None):
    target = Path(target).absolute()
    run_dir = Path(run_dir).absolute()
    index = run_dir / "miniprot.index"
    cmd = f"{miniprot} -t 16 -d {index} {target}"
    check = parallel_run_cmds(
        cmd, run_dir=run_dir, check="miniprot_index", check_path=check_path
    )
    return index



def miniprot_mapping(query, index, run_dir, check=False, onlycmd=False):
    query = Path(query).absolute()
    index = Path(index).absolute()
    run_dir = Path(run_dir).absolute()
    run = Path(run_dir)
    cmd = f"{miniprot} --gff-only -t 4 {index} {query} > {run}/miniprot.gff"
    if onlycmd:
        return cmd
    else:
        if check:
            check = ".miniprot_align"
            check = parallel_run_cmds(cmd, run_dir=run_dir, check=check)
        else:
            parallel_run_cmds(cmd, run_dir=run_dir)


def main():
    target = sys.argv[1]
    run_dir = sys.argv[2]
    query = sys.argv[3]
    index = miniprot_index(target, run_dir)
    miniprot_mapping(query, index, run_dir)


if __name__ == "__main__":
    main()
