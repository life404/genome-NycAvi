#! /usr/bin/python3

"""
Perform toga on each genes
"""
import pandas as pd
import re

from pathlib import Path
from concurrent import futures as cf
from modules.extract_miniprot_alignment import extract_region
from modules.lastz_on_miniprot import lastz_align
from modules.parallel_execute import run_cmd


class run_toga:
    def __init__(self, args):
        # Parsing the config file
        self.config = args.config

        # Obtaining the path of software
        self.toga = self._parse_config("Software", "toga")
        self.twoBitInfo = self._parse_config("Software", "ucsc_tools") / "twoBitInfo"

        # Obtaining the path of input file
        self.root = Path(args.workdir).resolve()
        self.target = self._parse_config("Path", "target")
        self.target_info = self._twoBitInfo(self.target)
        self.target_bed = self._parse_config("Path", "bed")
        self.target_isoform = self._parse_config("Path", "isoform")
        self.query = self._parse_config("Path", "query")
        self.query_info = self._twoBitInfo(self.query)
        self.ucsc_tools = self._parse_config("Software", "ucsc_tools")

        # Obtaining the value of additional parameters
        self.togaCBMermory = [
            int(i) for i in self._parse_config("Parameter", "togaCBMermory").split(",")
        ]
        self.togaParameter = self._parse_config("Parameter", "togaParameter")
        self.threads = args.threads

        self.bed_df = pd.read_csv(
            self.target_bed,
            sep="\t",
            header=None,
            names=[
                "chrom",
                "mRNAStart",
                "mRNAEnd",
                "name",
                "score",
                "strand",
                "exonStart",
                "exonEnd",
                "itemRgb",
                "cdsCount",
                "cdsSizes",
                "cdsStarts",
            ],
        )
        self.isoform_df = pd.read_csv(self.target_isoform, header=0, sep="\t")

        # Generating toga commands
        #self.cmds = self._generate_toga_cmd(gene_dir)

        # cmds_buff = open('./cmds', 'w')
        # cmds_buff.writelines(self.cmds)
        # cmds_buff.close()
        # self.cmds = [i for i in open('./cmds', 'r')]

    def _parse_config(self, first, second):
        return extract_region._parse_config(self, first, second)

    def _twoBitInfo(self, twoBit):
        return lastz_align._twoBitInfo(self, twoBit)

    def _prepre_toga_bed(self, gene, directory):
        iso_sub = self.isoform_df[self.isoform_df.GeneID == gene]
        iso_sub_path = directory / "toga.isoform"
        iso_sub.to_csv(iso_sub_path, sep="\t", header=False, index=False)

        bed_sub = self.bed_df[self.bed_df.name.isin(iso_sub.TransID.tolist())]
        bed_sub_path = directory / "toga.bed"
        bed_sub.to_csv(bed_sub_path, sep="\t", header=False, index=False)

        return iso_sub_path, bed_sub_path

    def _generate_toga_cmd(self, cmd_dir):
        cmds = list()
        print(f"Preparing the input files for toga......")
        for directory in cmd_dir:
            work_dir = Path(directory).resolve()
            toga_dir = work_dir.parent / "toga"
            toga_dir.mkdir(exist_ok=True, parents=True)
            nextflow_dir = toga_dir / "nextflow"

            """ Preparing the files for toga """
            gene = toga_dir.parts[-2]
            isoform_path, bed_path = self._prepre_toga_bed(gene, toga_dir)
            clean_chain = work_dir / "clean.chain"

            cmd = f"{self.toga} {clean_chain} {bed_path} {self.target} {self.query} --isoforms {isoform_path} --pd {toga_dir} --nd {nextflow_dir}  {self.togaParameter}"
            cmds.append(cmd)

        return cmds

    def check_toga(self, cmd_dir):
        """In sometime, the gene predicted by miniprot2 will not be regarded as PG in toga. So, in this study, the PG genes will be extracted and run this pipeline."""
        
        toga_results_df = pd.DataFrame(self._merge_toga_results(cmd_dir))
        toga_results_df = toga_results_df[(toga_results_df[2]=='PG') | (toga_results_df[2]=='E')]
        return toga_results_df[1].tolist()

    def _merge_toga_results(self, cmd_dir):            
        
        def parse_toga_results(directory):
            file_path = directory.parent / "toga" / "loss_summ_data.tsv"
            gene = directory.parts[-2]
            
            if file_path.is_file():
                return [i.strip().split('\t') for i in open(file_path, 'r') if i.startswith("GENE")]
            else:
                return [['GENE', gene, "E"]]
        
        print(f'Checkign the toga results......')
        merge_results = []
        with cf.ThreadPoolExecutor(max_workers=self.threads) as e:
            threads_lst = [e.submit(parse_toga_results, directory) for directory in cmd_dir]
            
            finish = 0
            all = len(cmd_dir)
            for threads in cf.as_completed(threads_lst):
                merge_results.extend(threads.result())
                print(f"{int(finish / all)}%", end='\r')
        print('Finish')
        return merge_results
    
    def results(self, gene_dir):
            work_dir = self.root
            check_point = work_dir / "checkpoint" / "toga.ok"
            continue_point = work_dir / "checkpoint" / "toga.continue"
            
            if continue_point.exists():
                gene_dir = [f"{i.strip()}/lastz" for i in open(continue_point, 'r')]

            minimum, maximum, step = self.togaCBMermory
            error_list = list()
            
            cmds = self._generate_toga_cmd(gene_dir)

            for max_memory in range(minimum, maximum, step):
                if cmds:
                    print(f"toga run on {max_memory}G memory")
                    cmds_tmp = [
                        f"{cmd} " + f"--cb {max_memory} --cesar_mem_limit {max_memory}"
                        for cmd in cmds
                    ]
                    shell = run_cmd(cmds_tmp, run_dir=work_dir)
                    results = shell.on_memory(
                        bar=True,
                        max_memory=max_memory,
                        threads=self.threads,
                        desc=f"toga: {max_memory}G",
                        unit="gene",
                    )

                    cmds = list()
                    rights = 0
                    for output, code, cmd in results:
                        if code != 0:
                            if (
                                "ALL genes require much more memory than the available"
                                in output
                            ):
                                cmds.append(re.sub(r"--cb.*", "", cmd))
                            else:
                                error_list.append(cmd)
                        else:
                            rights += 1

                    print(
                        f"There are {rights} gene run correctely, and {len(cmds)} genes need to be runed on more memory"
                    )
                else:
                    break

            print(
                f"Finally, there are ${len(error_list)} genes have error, please cheak them at {check_point}"
            )
            check_point.touch()
            check_buff = open(check_point, "w")
            check_buff.writelines([f"{i}\n" for i in error_list])
            check_buff.close()
