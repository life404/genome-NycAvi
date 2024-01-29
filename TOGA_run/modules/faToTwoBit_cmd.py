#! /usr/bin/python3

from pathlib import Path
from modules.parallel_execute import parallel_run_cmds
from modules.parallel_execute import shell_cmd_run
from configparser import ConfigParser
import sys

config = ConfigParser()
config.read(Path(__file__).absolute().parent / 'config.ini')
ucsc = config['Software']['ucsc_tools']
faToTwoBit_exe = Path(ucsc) / 'faToTwoBit'

def faToTwoBit(fa, run_dir, rmfa = False):
    
    fa = Path(fa).absolute()
    run_dir = Path(run_dir).absolute()
    twoBit = run_dir / f"{fa.name.replace('.fa', '.2bit')}"
    
    cmds = f"{faToTwoBit_exe} {fa} {twoBit}"
    output, code, command = shell_cmd_run(cmds, run_dir=run_dir)
    
    if rmfa:
        fa.unlink()
    return twoBit
        
def main():
    fa = Path(sys.argv[1])
    run_dir = Path(sys.argv[2])
    twoBit = faToTwoBit(fa, run_dir)
    print(twoBit)
    
if __name__ == "__main__":
    main()
	
        
        
    
    