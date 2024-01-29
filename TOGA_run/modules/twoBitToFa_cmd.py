#! /usr/bin/python3

from pathlib import Path
from modules.parallel_execute import parallel_run_cmds
from configparser import ConfigParser
import sys

config = ConfigParser()
config.read(Path(__file__).absolute().parent / "config.ini")
ucsc = config["Software"]["ucsc_tools"]
twoBitToFa_exe = Path(ucsc) / "twoBitToFa"


def twoBitToFa(twoBit, run_dir):
    twoBit = Path(twoBit)
    fa = run_dir / f"{twoBit.name.replace('.2bit', '.fa')}"
    check = f".twoBitToFa.{twoBit.name.replace('.2bit', '')}"
    cmds = f"{twoBitToFa_exe} {twoBit} {fa}"
    check = parallel_run_cmds(cmds, run_dir=run_dir, check=check)
    return check


def main():
    twoBit = Path(sys.argv[1])
    run_dir = Path(sys.argv[2])
    check = twoBitToFa(twoBit, run_dir)


if __name__ == "__main__":
    main()
