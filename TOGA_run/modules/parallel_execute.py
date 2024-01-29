#! /usr/bin/python3

import concurrent.futures as cf
import subprocess
import sys
import psutil
from queue import Queue
from tqdm import tqdm


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
            *args,
            **kwargs,
        )
        output, _ = run.communicate()
        code = run.wait()
        return output, code, cmd

    def on_parally(
        self,
        bar=False,
        desc="tqdm",
        unit="unit",
        threads=1,
        error="error.log",
        *args,
        **kwargs,
    ):
        threads = min(len(self.cmd_lst), threads) if len(self.cmd_lst) >= 1 else 1
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
                process_lst = [
                    e.submit(self._run, cmd, *args, **kwargs) for cmd in self.cmd_lst
                ]
                for process in cf.as_completed(process_lst):
                    output, code, cmd = process.result()
                    output_lst.append(output)
                    code_lst.append(code)
                    cmd_lst.append(cmd)
                    process_bar.update(1)
            process_bar.close()
        else:
            with cf.ProcessPoolExecutor(max_workers=threads) as e:
                process_lst = [
                    e.submit(self._run, cmd, *args, **kwargs) for cmd in self.cmd_lst
                ]
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

        return {"code": code_lst, "cmd": cmd_lst, "output": output_lst}

    def avilable_memory(self):
        mem = psutil.virtual_memory()
        available_memory_gb = mem.available / (1024**3)
        return int(available_memory_gb)

    def on_memory(
        self,
        bar=False,
        max_memory=5,
        free_memory=50,
        threads=64,
        desc="tqdm",
        unit="unit",
        *args,
        **kwargs,
    ):
        task_queue = Queue()
        for cmd in self.cmd_lst:
            task_queue.put((self._run, cmd, args, kwargs))

        task_lst = set()
        results_lst = set()
        pool = cf.ProcessPoolExecutor(max_workers=threads)

        if bar:
            
            process_bar = tqdm(
                total=len(self.cmd_lst), desc=f"{desc:25}", unit=f"{unit}"
            )
            while not task_queue.empty():
                if self.avilable_memory() > (int(free_memory) + int(max_memory)):
                    func, cmd, args, kwargs = task_queue.get()
                    task_lst.add(pool.submit(func, cmd, *args, **kwargs))

                
                done, task_lst = cf.wait(
                    task_lst, return_when="FIRST_COMPLETED", timeout=1
                )
                if done:
                    for process in done:
                        results_lst.add(process.result())
                        process_bar.update(1)
                else:
                    continue

            for process in cf.as_completed(task_lst):
                results_lst.add(process.result())
                process_bar.update(1)
            process_bar.close()

        else:
            while not task_queue.empty():
                if self.avilable_memory() > int(max_memory):
                    func, cmd, args, kwargs = task_queue.get()
                    task_lst.add(pool.submit(func, cmd, *args, **kwargs))

                done, task_lst = cf.wait(
                    task_lst, return_when="FIRST_COMPLETED", timeout=0.5
                )
                if done:
                    for process in done:
                        results_lst.add(process.result())
                else:
                    continue
            for process in cf.as_completed(task_lst):
                results_lst.add(process.result())
        return results_lst
