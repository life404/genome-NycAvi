#! /usr/bin/python3

from pathlib import Path
from configparser import ConfigParser
import sys

from modules.parallel_execute import parallel_run_cmds

config = ConfigParser()
config.read(Path(__file__).absolute().parent / "config.ini")
ucsc = config["Software"]["ucsc_tools"]
twoBitInfo_exe = Path(ucsc) / "twoBitInfo"


def twoBitInfo(twoBit, run_dir):
    run_dir = Path(run_dir).absolute()
    twoBit = Path(twoBit).absolute()
    info = run_dir / f"{twoBit.name.replace('.2bit', '.info')}"

    cmd = f"{twoBitInfo_exe} {twoBit} stdout | sort -k2rn > {info}"
    parallel_run_cmds(cmd, run_dir=run_dir)

    return info


def main():
    twoBit = sys.argv[1]
    run_dir = sys.argv[2]

    info = twoBitInfo(twoBit, run_dir)
    print(info)


if __name__ == "__main__":
    main()
