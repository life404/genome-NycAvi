#! /usr/bin/python3

from pathlib import Path
from configparser import ConfigParser
import shutil
import sys

from modules.parallel_execute import parallel_run_cmds
from modules.parallel_execute import shell_cmd_run

config = ConfigParser()
config.read(Path(__file__).absolute().parent / "config.ini")
ucsc = config["Software"]["ucsc_tools"]
lastz_exe = Path(ucsc) / "lastz"
axtChain_exe = Path(ucsc) / "axtChain"
RepeatFiller_exe = Path(config["Software"]["RepeatFiller"])
chainMergeSort_exe = Path(ucsc) / "chainMergeSort"


def lastz(
    target,
    query,
    run_dir,
    lastz_parameters="K=2400 H=2000 Y=9400 L=3000",
    error=None,
    error_path=None,
    cmdOnly=False,
):
    query = Path(query).absolute()
    target = Path(target).absolute()
    run_dir = Path(run_dir).absolute()
    axt = run_dir / f"{target.stem}-{query.stem}.axt"

    cmd = f"{lastz_exe} {target} {query} {lastz_parameters} --format=axt > {axt}"
    if cmdOnly:
        return cmd
    else:
        output, code, command = shell_cmd_run(cmd, run_dir=run_dir)
        if code != 0:
            if error_path and error:
                error_buff = open(error_path / error, 'w')
                error_buff.write(output[0])
                error.close()
            else:
                print(output[0])
        return axt


def axtChain(axt, target, query, run_dir, error=None, error_path=None, cmdOnly=False):
    axt = Path(axt).absolute()
    target = Path(target).absolute()
    query = Path(query).absolute()
    run_dir = Path(run_dir).absolute()
    chain = run_dir / f"{axt.stem}.chain"

    cmd = f"{axtChain_exe} -linearGap=loose {axt} {target} {query} {chain}"
    if cmdOnly:
        return cmd
    else:
        parallel_run_cmds(cmd, run_dir, error=error, error_path=error_path)
        output, code, command = shell_cmd_run(cmd, run_dir=run_dir)
        if code != 0:
            if error_path and error:
                error_buff = open(error_path / error, 'w')
                error_buff.write(output[0])
                error.close()
            else:
                print(output[0])
        return chain


def RepeatFiller(
    chain, target, query, run_dir, error=None, error_path=None, cmdOnly=False
):
    chain = Path(chain).absolute()
    target = Path(target).absolute()
    query = Path(query).absolute()
    run_dir = Path(run_dir).absolute()
    filled = run_dir / "filled.chain"

    cmd = f"{RepeatFiller_exe} -c {chain} -T2 {target} -Q2 {query} -v|grep -v '^$' > {filled}"
    if cmdOnly:
        return cmd
    else:
        output, code, command = shell_cmd_run(cmd, run_dir=run_dir)
        if "INFO:root:Found no new blocks to insert in this chain. Done!" in output[0] and code == 1:
            shutil.copy(chain, filled)
        #if code != 0:
        #    filled.touch()
        #    if error and error_path:
        #        error_buff = open(error_path / error, "w")
        #        error_buff.write("After RepeatFiller, the result is blank")
        #        error_buff.close()
        # if not output[0]:
        #    filled_chain_buff = open(filled, 'w')
        #    filled_chain_buff.writelines(
        #        [i.strip() for i in output[0] if i]
        #    )
        #    filled_chain_buff.close()
        # else:
        #    filled.touch()
        ##parallel_run_cmds(cmd, run_dir=run_dir)
        return filled
    
def chainMergeSort(chain_lst, run_dir, error = None, error_path = None):
    run_dir = Path(run_dir).absolute()
    
    chain_tmp = Path(run_dir) / 'tmp_chain.txt'
    chains_buff = open(chain_tmp, 'w')
    chains_buff.writelines([f"{i}\n" for i in chain_lst])
    chains_buff.close()

    cmd = f"{chainMergeSort_exe} -inputList={chain_tmp} > {run_dir}/target-query.chain"
    output, code, command = shell_cmd_run(cmd, run_dir=run_dir)
    if code != 0:
        if error and error_path:
            error_buff = open(error_path / error, 'w')
            error_buff.write(output[0])
            error_buff.close()
        else:
            print(output[0])
    
    chain_tmp.unlink()
    
    return run_dir / 'target-query.chain'

def main():
    target = sys.argv[1]
    query = sys.argv[2]
    run_dir = sys.argv[3]

    axt = lastz(target=target, query=query, run_dir=run_dir)
    print(axt)
    chain = axtChain(axt, target=target, query=query, run_dir=run_dir)
    print(chain)
    filled = RepeatFiller(chain=chain, target=target, query=query, run_dir=run_dir)
    print(filled)


if __name__ == "__main__":
    main()
