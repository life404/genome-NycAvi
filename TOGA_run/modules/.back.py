#! /usr/bin/python3

import concurrent.futures as cf
import subprocess
import sys
import psutil
from queue import Queue
from tqdm import tqdm


#def shell_cmd_run(command, run_dir):
#    run = subprocess.Popen(
#        command,
#        shell=True,
#        stdout=subprocess.PIPE,
#        stderr=subprocess.STDOUT,
#        text=True,
#        cwd=run_dir,
#    )
#    output, _ = run.communicate()
#    code = run.wait()
#    print(output)
#    return (output, code, command)
#
#
#def parallel_run_cmds(
#    cmds,
#    run_dir,
#    check=None,
#    check_path=None,
#    verbose=False,
#    num=16,
#    error=None,
#    error_path=None,
#):
#    if isinstance(cmds, str):
#        cmds = [cmds]
#
#    if check:
#        if check_path:
#            check_dir = Path(check_path)
#            check_dir.mkdir(exist_ok=True, parents=True)
#        else:
#            check_dir = Path(run_dir).absolute()
#
#        # check the check file to pass the run
#        check_file = check_dir / f"{check}.ok"
#
#        if check_file.exists():
#            print(f"Find the check point `{check_file}`, skip this step")
#        else:
#            num = min(len(cmds), num)
#
#            with cf.ProcessPoolExecutor(max_workers=num) as e:
#                process_lst = [
#                    e.submit(shell_cmd_run, command, run_dir) for command in cmds
#                ]
#            # Saving the stdout and error of each process
#            code_lst = list()
#            stdout_lst = list()
#            stderr_lst = list()
#
#            for process in cf.as_completed(process_lst):
#                output, code, command = process.result()
#                code_lst.append(code)
#
#                if code != 0:
#                    stderr_lst.append(f"The command:\n{command}\n")
#                    stderr_lst.append(output[0])
#                    stderr_lst.append(f"{'***'*30}\n")
#                stdout_lst.append(output[0])
#
#            if verbose:
#                for stdout in stdout_lst:
#                    for line in stdout:
#                        print(line)
#
#            if sum(code_lst) != 0:
#                if error:
#                    if not error_path:
#                        error_path = run_dir
#                    else:
#                        error_path = Path(error_path).absolute()
#                        print(
#                            f"The detailed error information can be checked in error log {error_path}/{error}.error"
#                        )
#                        error_buff = open(error_path / f"{error}.error", "w")
#                        error_buff.writelines(stderr_lst)
#                        error_buff.close()
#                else:
#                    print(
#                        f"Something is wrong in, Please set `verbose` true to check output"
#                    )
#            else:
#                check_file.touch()
#                print(f"Job(s) has finished, the check file is {check}")
#        return check_file
#
#    else:
#        with cf.ProcessPoolExecutor(max_workers=num) as e:
#            process_lst = [
#                e.submit(shell_cmd_run, command, run_dir) for command in cmds
#            ]
#
#        code_lst = list()
#        stdout_lst = list()
#        stderr_lst = list()
#
#        for process in cf.as_completed(process_lst):
#            output, code, command = process.result()
#            code_lst.append(code)
#
#            if code != 0:
#                stderr_lst.append(f"The command: \n{command}\n")
#                stderr_lst.append(output[0])
#                stderr_lst.append(f"{'***'*30}\n")
#            stdout_lst.append(output[0])
#
#        if verbose:
#            for stdout in stdout_lst:
#                for line in stdout:
#                    print(line)
#
#        if sum(code_lst) != 0:
#            if error:
#                if not error_path:
#                    error_path = run_dir
#                else:
#                    error_path = Path(error_path).absolute()
#                print(
#                    f"The detailed error information can be checked in error log {error_path}/{error}.error"
#                )
#                error_buff = open(error_path / f"{error}.error", "w")
#                error_buff.writelines(stderr_lst)
#                error_buff.close()
#            else:
#                print(
#                    f"Something is wrong in, Please set `verbose` true to check output"
#                )
#

class run_cmd:
    """run the shell commands"""

    def __init__(self, cmd, run_dir):
        if isinstance(cmd, str):
            self.cmd_lst = [cmd]
        else:
            self.cmd_lst = cmd

        self.run_dir = run_dir

    def _run(self, cmd, *args, **kwargs):
        run = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=self.run_dir,
            *args, **kwargs 
        )
        output, _ = run.communicate()
        code = run.wait()
        return output, code, cmd

    def on_parally(self, bar=False, desc="tqdm", unit="unit", threads=1, error = 'error.log', out = False, *args, **kwargs):
        threads = min(len(self.cmd_lst), threads)
        # collecting the output
        code_lst = list()
        output_lst = list()
        cmd_lst = list()

        # ask wheather show the process bar
        if bar:
            process_bar = tqdm(
                total=len(self.cmd_lst), desc=f"{desc:25}", unit=f"{unit}"
            )
            with cf.ProcessPoolExecutor(max_workers=threads) as e:
                process_lst = [e.submit(self._run, cmd, *args, **kwargs) for cmd in self.cmd_lst]
                for process in cf.as_completed(process_lst):
                    output, code, cmd = process.result()
                    output_lst.append(output)
                    code_lst.append(code)
                    cmd_lst.append(cmd)
                    process_bar.update(1)
            process_bar.close()
        else:
            with cf.ProcessPoolExecutor(max_workers=threads) as e:
                process_lst = [e.submit(self._run, cmd, *args, **kwargs) for cmd in self.cmd_lst]
                for process in cf.as_completed(process_lst):
                    output, code, cmd = process.result()
                    output_lst.append(output)
                    code_lst.append(code)
                    cmd_lst.append(cmd)

        # ask wheather have a error
        if sum(code_lst) != 0:
            print(f"Something is wrong, please check the error log")
            error_log = self.run_dir / error
            error_buff = open(error_log, "w")
            error_buff.writelines(output_lst)
            error_buff.close()
        
        if out:
            return {'code':code_lst, 'cmd':cmd_lst, 'output':output_lst}
    
    def avilabel_memory(self):
        mem = psutil.virtual_memory()
        available_memory_gb = mem.available / (1024**3)
        return available_memory_gb - 20
    
    def process_pool(self, max_memory):
        available_memory_gb = self.avilabel_memory()
        num = int(available_memory_gb / max_memory)
        pool = cf.ProcessPoolExecutor(max_workers=num)
        return pool
        
    
        
    #def on_memory(self, bar=False, max_memory=5,desc='tqdm', unit='unit', error='error.log', out=False, *args, **kwargs):
    #    
    #    tasks = Queue()
    #    for cmd in self.cmd_lst:
    #        tasks.put(
    #            self._run,
    #            cmd = cmd,
    #            *args,
    #            **kwargs
    #        )
    #        
    #    if bar:
    #        process_bar = tqdm(total = len(self.cmd_lst), desc=f"{desc:25}", unit=f"{unit}")
    #        
    #        pool = self.process_pool(max_memory)
    #        
    #        while not tasks.empty():
    #            if self.avilabel_memory > (max_memory * )
        
