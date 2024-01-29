#! /usr/bin/python3

import psutil
import shutil
from tqdm import tqdm
from pathlib import Path
from concurrent import futures as cf
from queue import Queue

from modules.lastz_cmd import lastz, axtChain, RepeatFiller, chainMergeSort
from modules.chainCleaner_cmd import chainCleaner
from modules.query_split import query_split_by_chromosome


class CustomExecutor_on_cpu_load:
    """
    定义一个类，用于检测 CPU 的负载，基于 CPU 的负载动态调整任务的提交
    每次提交任务前，检测当前 CPU 的负载，如果当前空闲的 CPU大于阈值，那么提交任务，然后再次判断 CPU 的负载，如果空闲小于阈值，那么不在提交任务
    由于我们将query 物种的基因组基于染色体进行了拆分，进行并行的运算，因此对于每一个基因来说，阈值就是每个 query 物种基因组的序列的数目。但是如果某些染色体的长度非常短，那么 lastz 运行也不会消耗太多时间。因此，我们规定大于 128k 2bit的文件进行统计
    """

    def __init__(self, cpu_threshold):
        self.pool = cf.ProcessPoolExecutor()
        self.tasks = Queue()
        self.results = list()
        self.cpu_threshold = cpu_threshold

    def submit(self, func, *args, **kwargs):
        # 构建 queue 队列，存储任务的参数
        self.tasks.put((func, args, kwargs))

    def idle_cpu(self):
        # cpu 判断
        cpu_percent_lst = psutil.cpu_percent(interval=1, percpu=True)
        idle_cpu = [cpu for cpu in cpu_percent_lst if cpu < 20]
        return len(idle_cpu)

    def execute(self):
        process_bar = tqdm(
            total=self.tasks.qsize(), desc="lastz align again", unit="gene"
        )
        tasks_lst = set()
        while not self.tasks.empty():
            # cpu 检测，动态提交任务
            if self.idle_cpu() > self.cpu_threshold:
                func, args, kwargs = self.tasks.get()
                tasks_lst.add(self.pool.submit(func, *args, **kwargs))
            else:
                pass
            # 检测任务是否完成，更新进度条
            done, tasks_lst = cf.wait(
                tasks_lst, return_when="FIRST_COMPLETED", timeout=1
            )
            if done:
                for process in done:
                    self.results.append(process.result())
                    process_bar.update(1)
            else:
                continue

        for process in cf.as_completed(tasks_lst):
            self.results.append(process.result())
            process_bar.update()
        process_bar.close()

        return self.results


def lastz_rerun_pipeline1(target, query, run_dir):
    run_dir = Path(run_dir).absolute()
    if not run_dir.exists():
        run_dir.mkdir(exist_ok=True, parents=True)

    axt = lastz(target, query, run_dir=run_dir)
    chain = axtChain(axt=axt, target=target, query=query, run_dir=run_dir)
    return chain


def lastz_rerun_pipeline2(chain_lst, target, query, run_dir, error, error_path):
    run_dir = Path(run_dir).absolute()
    if not run_dir.exists():
        run_dir.mkdir(exist_ok=True, parents=True)

    chain = chainMergeSort(
        chain_lst=chain_lst,
        run_dir=run_dir,
        error=f"chainMergeSort.{error}",
        error_path=error_path,
    )
    filled = RepeatFiller(
        chain=chain,
        target=target,
        query=query,
        run_dir=run_dir,
        error=f"RepeatFiller.{error}",
        error_path=error_path,
    )
    clean = chainCleaner(
        chain=filled,
        target=target,
        query=query,
        run_dir=run_dir,
        error=f"RepeatFiller.{error}",
        error_path=error_path,
    )
    return clean


def lastz_rerun(gene, query_2bit, query_split_2bit_lst, root_dir, error_path=None):
    run_dir = Path(root_dir).absolute() / f"{gene}" / "lastz"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    
    query_2bit = Path(query_2bit).absolute()
    axt_dir = run_dir / "axt"
    target_2bit = run_dir.parent / "target.2bit"
    
    # replace the old 2bit file
    query_2bit_raw = run_dir.parent / "query.2bit"
    query_2bit_raw.unlink(missing_ok=True)
    query_2bit_raw.symlink_to(query_2bit)

    error_path = Path(error_path).absolute()

    with cf.ProcessPoolExecutor() as e:
        process_lst = [
            e.submit(
                lastz_rerun_pipeline1,
                target=target_2bit,
                query=query_split_2bit,
                run_dir=axt_dir,
            )
            for query_split_2bit in query_split_2bit_lst
        ]
        chain_lst = list()
        for process in cf.as_completed(process_lst):
            chain = process.result()
            chain_lst.append(chain)

    clean = lastz_rerun_pipeline2(
        chain_lst=chain_lst,
        target=target_2bit,
        query=query_2bit,
        run_dir=run_dir,
        error=gene,
        error_path=error_path,
    )

    return clean.parents[1]

def lastz_rerun_on_cpu_percent(gene_lst, query, root_dir, error_path=None):
    query = Path(query).absolute()
    query_split = query_split_by_chromosome(query)
    query_2bit = query_split / f"{query.stem}.2bit"
    query_split_2bit_lst = [
        file
        for file in list(query_split.glob("*.2bit"))
        if file.stem != f"{query.stem}"
    ]

    error_path = Path(error_path).absolute()
    cpu_threshold = len(
        [i for i in query_split_2bit_lst if i.stat().st_size > 1024 * 128]
    )
    custom_e = CustomExecutor_on_cpu_load(cpu_threshold=cpu_threshold)
    for gene in gene_lst:
        custom_e.submit(
            lastz_rerun,
            gene=gene,
            query_2bit=query_2bit,
            query_split_2bit_lst=query_split_2bit_lst,
            root_dir=root_dir,
            error_path=error_path,
        )
    results_lst = custom_e.execute()

    return results_lst
