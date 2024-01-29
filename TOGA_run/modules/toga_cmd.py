#! /usr/bin/python3

from pathlib import Path
from configparser import ConfigParser
import sys
import concurrent.futures as cf
from tqdm import tqdm
import psutil
import time
from queue import Queue

from modules.parallel_execute import shell_cmd_run

config = ConfigParser()
config.read(Path(__file__).absolute().parent / "config.ini")
toga_exe = Path(config["Software"]["toga"])

class CustomeExecutor_on_avilabel_memory:
    
    def __init__(self, cmx, cjn):
        self.tasks = Queue()
        self.results = list()
        self.cmx = int(cmx)
        self.cjn = int(cjn)
    
    def avilabel_memory(self):
        mem = psutil.virtual_memory()
        available_memory_gb = mem.available / (1024**3)
        return (available_memory_gb - 100)
    
    def process_pool(self):
        mem = psutil.virtual_memory()
        available_memory_gb = mem.available / (1024**3)
        num = int(available_memory_gb / self.cmx)
        pool = cf.ProcessPoolExecutor(max_workers=num)
        return pool
    
    def submit(self, func, *args, **kwargs):
        self.tasks.put((func, args, kwargs))
       
    def execute(self):
        process_bar = tqdm(total = self.tasks.qsize(), desc=f"{'TOGA':25}", unit='gene')
        tasks_lst = set()
        pool = self.process_pool()
        
        while not self.tasks.empty():
            if self.avilabel_memory() > (self.cmx * self.cjn):
                func, args, kwargs = self.tasks.get()
                tasks_lst.add(pool.submit(func, *args, **kwargs))
            
            done, tasks_lst = cf.wait(tasks_lst, return_when='FIRST_COMPLETED', timeout=1)
            if done:
                for process in done:
                    self.results.append(process.result())
                    process_bar.update(1)
            else:
                continue
        
        for process in cf.as_completed(tasks_lst):
            self.results.append(process.result())
            process_bar.update(1)
        process_bar.close()
        
        return self.results
                

def toga(run_dir, error_path=None, chn=1, cjn=1, cmx=0):
    run_dir = Path(run_dir).absolute()
    chain = run_dir / "lastz" / "clean.chain"
    bed = run_dir / "toga.bed"
    query = run_dir / "query.2bit"
    target = run_dir / "target.2bit"
    isoform = run_dir / "isoform.tsv"
    if str(run_dir).endswith("rerun"):
        error = f"toga.rerun.{run_dir.parts[-2]}"
    else:
        error = f"toga.{run_dir.stem}"

    toga_dir = run_dir / "toga"
    toga_dir.mkdir(exist_ok=True, parents=True)

    cmd = f"{toga_exe} {chain} {bed} {target} {query} -i {isoform} --kt --pd {toga_dir} --nd {toga_dir}/nextflow --chn {chn} --cjn {cjn} --ms --cb {cmx} --cesar_mem_limit {cmx}"
    output, code, command = shell_cmd_run(cmd, run_dir=run_dir)
    # test_buff = open(error_path / f"{error}.test", "w")
    # test_buff.writelines([f"{code}\n", output[0]])
    # test_buff.close()

    if code == 0:
        flag = "ok"
    else:
        if "ALL genes require much more memory than the available" in output[0]:
            flag = "more_mem"
        else:
            flag = "error"
            error_buff = open(error_path / f"{error}.error", "w")
            error_buff.writelines(output[0])
            error_buff.close()
    return flag, run_dir


def toga_base_on_memory(run_dir_lst, error_path, cmx, cjn):
    toga_finish_lst = list()
    toga_rerun_lst = list()
    toga_error_lst = list()
    
    custome_e = CustomeExecutor_on_avilabel_memory(cmx=cmx, cjn=cjn)
    for run_dir in run_dir_lst:
        custome_e.submit(
            toga,
            run_dir = run_dir,
            error_path = error_path,
            cmx = cmx,
            cjn = cjn
        )
    results = custome_e.execute()
    for result in results:
        if result[0] == 'ok':
            toga_finish_lst.append(Path(result[1]) / "toga" / "loss_summ_data.tsv")
        elif result[0] == 'more_mem':
            toga_rerun_lst.append(result[1])
        else:
            toga_error_lst.append(result[1])
    return toga_finish_lst, toga_rerun_lst, toga_error_lst

    #def check_available_memory(cmx, cjn):
    #    mem = psutil.virtual_memory()
    #    availabel_memory_gb = mem.available / (1024**3)
    #    return availabel_memory_gb >= 100 + cmx * cjn

    #def threads_num(cmx):
    #    mem = psutil.virtual_memory()
    #    available_memory_gb = mem.available / (1024**3)
    #    return int(available_memory_gb / cmx)

    #threads = threads_num(cmx)
    #with cf.ProcessPoolExecutor(max_workers=threads) as e:
    #    process_lst = list()
    #    for run_dir in run_dir_lst:
    #        while not check_available_memory(cmx, cjn):
    #            time.sleep(0.5)
    #        process_lst.append(e.submit(toga, run_dir, error_path, cmx=cmx, cjn=cjn))
    #        

    #    for process in cf.as_completed(process_lst):
    #        flag, dir = process.result()
    #        if flag == "ok":
    #            toga_finish_lst.append(Path(dir) / "toga" / "loss_summ_data.tsv")
    #        elif flag == "more_mem":
    #            toga_rerun_lst.append(dir)
    #        else:
    #            toga_error_lst.append(dir)
    #        process_bar.update(1)
    #return toga_finish_lst, toga_rerun_lst, toga_error_lst


def main():
    run_dir = sys.argv[1]

    toga_dir = toga(run_dir)
    print(toga_dir)


if __name__ == "__main__":
    main()
    re.match()