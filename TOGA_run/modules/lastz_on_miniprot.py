#! /usr/bin/python3

"""
In this step, the region infomation from previous step will be used to run the lastz aling. The pipeline of lastz align is referenced on the make_lastz_chains. If a coding gene of target can be aliged to query genome sucessfully, the regions of coding gene from target and potential region from query genome will be used as input in alignemtn of lastz. Miniprot may fail to detect coding regions of certain genes in the queried genome due to its limited detection capability or the absence of the gene. For these genes that not be detected, we will align these genes to the whole query genome. 
"""

from pathlib import Path
import random
import string
import pandas as pd
from modules.parallel_execute import run_cmd
from modules.extract_miniprot_alignment import extract_region
from concurrent import futures as cf
import os
import shutil
import argparse
from tqdm import tqdm
from datetime import datetime

class lastz_align:
    """The class defined the pipeline of lastz alignment"""

    def __init__(self, args):
        # Inheriting config parsing function from extract_region
        # self._parse_config = extract_region._parse_config()

        """Parsing the config file and obtain values"""

        # values of software
        self.config = args.config
        self.lastz = self._parse_config("Software", "ucsc_tools") / "lastz"
        self.axtChain = self._parse_config("Software", "ucsc_tools") / "axtChain"
        self.chainMergeSort = (
            self._parse_config("Software", "ucsc_tools") / "chainMergeSort"
        )
        self.chainSort = self._parse_config("Software", "ucsc_tools") / "chainSort"
        self.chainScore = self._parse_config("Software", "ucsc_tools") / "chainScore"
        self.RepeatFiller = (
            self._parse_config("Software", "ucsc_tools") / "RepeatFiller.py"
        )
        self.chainExtractID = (
            self._parse_config("Software", "ucsc_tools") / "chainExtractID"
        )
        self.chainCleaner = (
            self._parse_config("Software", "ucsc_tools") / "chainCleaner"
        )
        self.chainNet = self._parse_config("Software", "ucsc_tools") / "chainNet"
        self.twoBitInfo = self._parse_config("Software", "ucsc_tools") / "twoBitInfo"

        # parameters of softwares

        self.lastz_paramerter = " ".join(
            [
                f'H={self._parse_config("Parameter", "lastz_h")}',
                f'Y={self._parse_config("Parameter", "lastz_y")}',
                f'L={self._parse_config("Parameter", "lastz_l")}',
                f'K={self._parse_config("Parameter", "lastz_k")}',
            ]
        )
        self.chainMinScore = self._parse_config("Parameter", "chainMinScore")
        self.chainLinearGap = self._parse_config("Parameter", "chainLinearGap")
        self.cleanChainParameter = self._parse_config(
            "Parameter", "cleanChainParameter"
        )
        self.repeatFillerParameter = self._parse_config(
            "Parameter", "repeatFillerParamerter"
        )

        # values of path
        self.root = Path(args.workdir).resolve()
        self.root.mkdir(exist_ok=True, parents=True)
        self.check = self.root / 'checkpoint'
        self.check.mkdir(exist_ok=True, parents=True)
        self.log = self.root / 'log'
        self.root.mkdir(exist_ok=True, parents=True)
        
        self.target = self._parse_config("Path", "target")
        self.target_info = self._twoBitInfo(self.target)
        self.query = self._parse_config("Path", "query")
        self.query_info = self._twoBitInfo(self.query)
        self.ucsc_tools = self._parse_config("Software", "ucsc_tools")

        # other
        self.threads = args.threads

        # Generate commands of lastz pipeline and saving them in a bash file """
        # In this programe, we put the genes that align to whole query in the front of genes that align to specific region of query. Beacuse, genes that align to whole query will cost more time that align to specific region of query.

    def _parse_config(self, first, second):
        return extract_region._parse_config(self, first, second)

    def __shuffle_name(self, k):
        name_lst = list()
        while len(set(name_lst)) < k:
            name_lst.append(
                "".join(random.sample(string.ascii_lowercase + string.digits, 4))
            )
        return list(set(name_lst))

    """
    重写生成lastz命令, 合并两种lastz命令生成方式，完全基于miniprot的结果进行命令的生成。对于miniprot无法比对的区域，基于上下文环境进行了判断， 不再使用全基因组比对，节省比对时间
    
    为每个基因分配4个线程进行运算
    如果单个比对区域，无需使用chainMergeSort, 因为所有lastz运算生成随机名字，因此即使使用单个比对区域，也使用chainMergeSort生成最终文件，这么做并不会明显增加运行成本
    
    """

    def _lastz(self, directory, region):
        work_dir = Path(directory)
        
        # 防止lastz报错，将所有位置信息都转换为整数
        target_regions = region.apply(
            lambda x: f"\"{x['chr_bed']}[{int(x['start_bed'])}, {int(x['end_bed'])}][multiple]\"",
            axis=1,
        ).tolist()

        # 判断，是否存在比对区域，如果不存在比对区域，那么和整条染色体进行比对
        query_regions = region.apply(
            lambda x: f"\"{x['chr_miniprot']}[{int(x['start_miniprot'])}, {int(x['end_miniprot'])}][multiple]\""
            if not pd.isna(x["start_miniprot"])
            else f"\"{x['chr_miniprot']}[multiple]\"",
            axis=1,
        ).tolist()

        shuffle_names = self.__shuffle_name(len(query_regions))
        cmd = f"mkdir -p {work_dir}/tmp_lastz\n"

        # 即使設置parallel線程爲4，如果命令只有一條也是單線程運算
        cmd += f'parallel -j 4 --xapply "lastz "{self.target}/{{1}}" "{self.query}/{{2}}" {self.lastz_paramerter} --format=axt | axtChain -linearGap=loose -minScore=1000 stdin {self.target} {self.query} {work_dir}/tmp_lastz/{{3}}.chain" ::: {" ".join(target_regions)} ::: {" ".join(query_regions)} ::: {" ".join(shuffle_names)}\n'
        cmd += self._chainMergeSort(work_dir)
        cmd += f"rm -rf {work_dir}/tmp_lastz\n"

        return cmd

    def _axtChain(self, directory):
        """Generating axtChain commands"""
        work_dir = Path(directory)
        cmd = f"{self.axtChain} -minScore={self.chainMinScore} -linearGap={self.chainLinearGap} {work_dir}/lastz.axt {self.target} {self.query} {work_dir}/lastz.chain"
        return cmd

    def _chainMergeSort(self, directory):
        cmd = f"find {directory}/tmp_lastz -name '*.chain' | chainMergeSort -inputList=stdin > {directory}/lastz.chain\n"
        return cmd

    def _fillChain(self, directory):
        """Generating RepeatFill.py commands"""
        work_dir = Path(directory)
        cmd = f"mkdir -p {work_dir}/repeatfill_temp"
        cmd += f"\n{self.RepeatFiller} --workdir {work_dir}/repeatfill_temp --chainExtractID {self.chainExtractID} --lastz {self.lastz} --axtChain {self.axtChain} --chainSort {self.chainSort} -c {work_dir}/lastz.chain -T2 {self.target} -Q2 {self.query} {self.repeatFillerParameter} | {self.chainScore} -linearGap={self.chainLinearGap} stdin {self.target} {self.query} stdout | {self.chainSort} stdin {work_dir}/fill.chain\n"
        cmd += f"rm -rf {work_dir}/repeatfill_temp\n"
        return cmd

    def _chainCleaner(self, directory):
        """Generating chainCleaner command"""
        work_dir = Path(directory)
        cmd = f"{self.chainCleaner} {work_dir}/fill.chain {self.target} {self.query} {work_dir}/clean.chain {work_dir}/removedSuspects.bed -linearGap={self.chainLinearGap} -tSizes={self.target_info} -qSizes={self.query_info} {self.cleanChainParameter}\n"
        return cmd

    def _twoBitInfo(self, twoBit):
        twoBit = Path(twoBit).resolve()
        info = twoBit.parent / f"{twoBit.stem}.info"
        if not info.exists():
            cmd = f"{self.twoBitInfo} {twoBit} stdout"
            execute = run_cmd(cmd, run_dir=self.workdir)
            out = execute.on_parally()
        else:
            pass
        return info

    """
    不再使用两个生成命令的方式，进行合并，仍然區分多線程任務，
    """

    def _generate_cmds(self, region_info):

        # 对region_info，进行分割，如果一个基因只有一个比对的区域，那么优先运行，保证运行效率。如果一个基因比对多个区域，那么分配四个线程，总线程数除以4，防止超出最大线程
        job_s, job_m = list(), list()
        for _, region in region_info.groupby("GeneID", as_index=False):
            cmd_lst = list()
            directory = self.root / "work" / region["GeneID"].iloc[0] / "lastz"
            directory.mkdir(exist_ok=True, parents=True)
            
            cmd_lst.append(self._lastz(directory, region))
            cmd_lst.append(self._fillChain(directory))
            cmd_lst.append(self._chainCleaner(directory))
            bash_buff = open(directory / "run_lastz.sh", "w")
            bash_buff.writelines(cmd_lst)
            bash_buff.close()
            
            if len(region) == 1:
                job_s.append(directory)
            else:
                job_m.append(directory)
                            
        return job_s, job_m

    def _perform_cmds(self, cmd_dir, threads=None):
        """Setting the environment variable"""
        env = os.environ.copy()
        env["PATH"] += f":{self.ucsc_tools}"

        if threads is None:
            threads = self.threads

        """Running the bash file"""
        cmds = [f"bash {i}/run_lastz.sh" for i in cmd_dir]
        shell = run_cmd(cmds, run_dir=self.root / "work")
        out = shell.on_parally(
            bar=True,
            threads=threads,
            desc="lastz",
            unit="gene",
            env=env,
        )

        return cmd_dir

    def _return_dirs(self):
        return self.cmd_dirs

    def _blank_chain_file(self, chain):
        chain_buff = [i for i in open(chain, "r") if not i.startswith("#")]
        return chain if not chain_buff else None

    def run(self, region_info):
        
        # run cmd,不同的任务使用不同的线程
        print(f"Generating the lastz commands ......")
        start_time = datetime.now()
        job_s, job_m = self._generate_cmds(region_info)
        end_time = datetime.now()
        print(f"{'time used':25}: {end_time - start_time}")

        print('Runing lastz alignment ......')
        start_time = datetime.now()
        check_path = self.check / "lastz.ok"
        if check_path.exists():
            cmd_dir = job_s + job_m
        else:
            cmd_dir_s = self._perform_cmds(job_s, self.threads)
            cmd_dir_m = self._perform_cmds(job_m, int(self.threads / 4))
            cmd_dir = cmd_dir_s + cmd_dir_m
            check_path.touch()
        end_time = datetime.now()
        print(f"{'time used':25}: {end_time - start_time}")

        # check the lastz chain. Somethimes, the region extracted from miniprot may produce a blank lastz results. These genes with blank chain files will be detected and re-run alignment against to whole query genome again
        # 理论上不应该，所有的gene应该都可以比对上，因为根据上下文环境找了潜在的比对范围，但是不排除特殊情况
        print("Final checking......")
        with cf.ThreadPoolExecutor() as e:
            process_lst = [
                e.submit(self._blank_chain_file, chain)
                for chain in [i / "lastz.chain" for i in cmd_dir]
            ]
            finish = 0
            all = len(cmd_dir)
            blank_chain_dir = list()
            for process in cf.as_completed(process_lst):
                finish += 1
                if process.result() != None:
                    blank_chain_dir.append(process.result())
                percent = int(finish / all * 100)
                print(f"{percent}%", end="\r")
        blank_chain_gene = [f"{Path(i).parts[-3]}" for i in blank_chain_dir]
        print(
            f"Finish checking, there are {len(blank_chain_gene)} genes with blank lastz output (lastz.chain)"
        )
        print(*blank_chain_gene, sep='\n')
#        blank_region_info = self.region_info.copy()
#        blank_region_info = blank_region_info[
#            blank_region_info.GeneID.isin(blank_chain_gene)
#        ]
#        rerun_cmd_dir = self._generate_cmds2(blank_region_info)
#        self._perform_cmds(rerun_cmd_dir, int(self.threads / 4))
#
#        # Sometimes, the RepeatFiller will not print any results. We will used the raw chain file of lastz to replace the chain of RepeatFiller and chainCleaner
        bad_chain = [
            Path(i) for i in cmd_dir if Path(i / "clean.chain").stat().st_size == 0
        ]
        print(
            f"There are {len(bad_chain)} genes can't pass the RepeatFiller, which will be replace with the raw lastz.chain. The informat have been saved in {check_path}"
        )
        # 保存bad chain 至log文件
        for i in bad_chain:
            raw = Path(i) / "lastz.chain"
            clean = Path(i) / "clean.chain"
            shutil.copy(raw, clean)
        error_log = self.log / 'lastz.log'    
        check_buff = open(error_log, "w")
        check_buff.writelines([f"{i.parts[-2]}\t{i / 'clean.chain'}\n" for i in bad_chain])
        check_buff.close()
        
        return cmd_dir


def make_args():
    parse = argparse.ArgumentParser()
    parse.add_argument(
        "-c", "--config", type=str, dest="config", help="The config file"
    )
    parse.add_argument(
        "-t",
        "--threads",
        type=int,
        dest="threads",
        default=32,
        help="How many processes were used in parally run. default[32]",
    )
    parse.add_argument(
        "-o",
        "--output",
        type=str,
        dest="workdir",
        help="The directory was used to save the results",
    )
    args = parse.parse_args()
    return args


def main():
    args = make_args()
    # step_one = extract_region(args)
    # region_info = step_one.results()
    # region_info.to_csv('/home/panda2bat/Avivorous_bat/script/TOGA_run/modules/region_info.csv', sep = '\t', index = False, header = True)
    region_info = pd.read_csv(
        "/home/panda2bat/Avivorous_bat/script/TOGA_run/test/output.new/lastz_region",
        sep="\t",
    )
    step_two = lastz_align(args, region_info.sample(n=500))


if __name__ == "__main__":
    main()
