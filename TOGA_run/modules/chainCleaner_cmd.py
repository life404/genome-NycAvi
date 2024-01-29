#! /usr/bin/python3

from pathlib import Path
from configparser import ConfigParser
import sys
import subprocess
import os

from modules.twoBitInfo_cmd import twoBitInfo

config = ConfigParser()
config.read(Path(__file__).absolute().parent / "config.ini")
ucsc = config["Software"]["ucsc_tools"]
chainCleaner_exe = Path(ucsc) / "chainCleaner"

# set the environment variable
env = os.environ.copy()
env["PATH"] += f":{ucsc}"
env["PATH"] += f":/home/panda2bat/TOOLS/GenomeAlignmentTools/src"


def shell_run_cmd(command, run_dir, env=None):
    run = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        cwd=run_dir,
        env=env,
    )
    output = run.communicate()
    code = run.wait()
    return (output, code, command)


def chainCleaner(chain, target, query, run_dir, error=None, error_path=None):
    chain = Path(chain).absolute()
    target = Path(target).absolute()
    query = Path(query).absolute()
    run_dir = Path(run_dir).absolute()

    targetInfo = twoBitInfo(target, run_dir=run_dir)
    queryInfo = twoBitInfo(query, run_dir=run_dir)

    clean = run_dir / "clean.chain"

    cmd = f"{chainCleaner_exe} {chain} {target} {query} {clean} {run_dir}/removeSuspects.bed -tSizes={targetInfo} -qSizes={queryInfo} -doPairs -minBrokenChainScore=75000 -linearGap=loose"
    output, code, command = shell_run_cmd(cmd, run_dir=run_dir, env=env)

    if code != 0:
        print(output)
        clean.touch()
        if error and error_path:
            error_buff = open(error_path / error, "w")
            error_buff.write("The cleaned chain is blank")
            error_buff.close()
    else:
        return clean


def main():
    chain = sys.argv[1]
    target = sys.argv[2]
    query = sys.argv[3]
    run_dir = sys.argv[4]

    clean = chainCleaner(chain=chain, target=target, query=query, run_dir=run_dir)
    print(clean)


if __name__ == "__main__":
    main()
