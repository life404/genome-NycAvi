#! /usr/bin/python3

from pathlib import Path
from modules.parallel_execute import parallel_run_cmds
from modules.parallel_execute import shell_cmd_run
from configparser import ConfigParser
import sys

config = ConfigParser()
config.read(Path(__file__).absolute().parent / "config.ini")
ucsc = config["Software"]["ucsc_tools"]
chainMergeSort_exe = Path(ucsc) / "chainMergeSort"

def chainMergeSort(inputList, run_dir, error=None, error_path=None):
    inputList = Path(inputList).absolute()
    run_dir = Path(run_dir).absolute()

    cmd = f"{chainMergeSort_exe} {inputList} > {run_dir}/target-all_query.chain"
    output, code, command = shell_cmd_run(cmd, run_dir)

    return run_dir / "target-all_query.chain"
