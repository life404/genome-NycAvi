#! /usr/bin/python3
"""
The aim of this step is to identify the potential region of gene coding region. Miniprot can detect the coding gene from genome based on the orthologous protein fastly. Therefore, we used miniprot to detect the gene coding region from genome. Second, Next, we extracted the coding region of the gene, as well as the flanking regions 10mb upstream and downstream. These regions were be used in lastz alignemnt. 
"""

from pathlib import Path
from Bio import SeqIO
from configparser import ConfigParser
import argparse
import numpy as np
import pandas as pd
import math
from concurrent import futures as cf
from tqdm import tqdm
import re
import time

# manual modules
from modules.parallel_execute import run_cmd


class extract_region:
    """Obtainnin the potentail lastz region based on the gene-level"""

    def __init__(self, args):
        # parsing the config file
        self.config = Path(args.config)

        # getting the path of miniprot
        self.miniprot = self._parse_config("Software", "miniprot")
        self.twoBitToFa = self._parse_config("Software", "ucsc_tools") / "twoBitToFa"
        self.twoBitInfo = self._parse_config("Software", "ucsc_tools") / "twoBitInfo"

        # Making the root directory of results
        self.root = Path(args.work_dir).resolve()
        self.root.mkdir(exist_ok=True, parents=True)
        # Making the check directory, under the root directory
        self.check = self.root / "checkpoint"
        self.check.mkdir(exist_ok=True, parents=True)
        # Making the work directory, under the root directory
        self.work_dir = self.root / "work"
        self.work_dir.mkdir(exist_ok=True, parents=True)

        self.target = self._parse_config("Path", "target")
        self.target_bed = self._parse_config("Path", "bed")
        self.target_faa = self._parse_config("Path", "faa")
        self.target_gff = self._parse_config("Path", "gff")
        self.target_isoform = self._parse_config("Path", "isoform")

        self.query = self._parse_config("Path", "query")
        self.query_f = self._transform_2bit_to_fa(
            self.query, Path(f"{self.query.parent}/{self.query.stem}.fa")
        )

        self.select = self._parse_config("Path", "select")
        self.qplank = self._parse_config("Parameter", "qplank")
        self.tplank = self._parse_config("Parameter", "tplank")

        # other
        self.threads = args.threads

        """Perfomming miniprot analysis"""

        self.miniprot_index = self._miniprot_index()

        # Parsing the target files and merge into one df
        self.target_df = self._read_target()

        # run the miniprot alignment
        self.target_df = self._miniprot_alignment()

        # Base on the alignment split the query file
        self.lastz_align_region = self._extract_lastz_align_region()

    def _parse_config(self, first, second):
        def ispath(x):
            return Path(x).is_file() or Path(x).is_dir()

        config_parse = ConfigParser()
        config_parse.optionxform = lambda option: option
        config_parse.read(self.config)
        value = config_parse[first][second]
        if value.isnumeric():
            return int(value)
        elif ispath(value):
            return Path(value)
        else:
            return str(value)

    def _read_target(self):
        """This function was used to obtain the max intron length"""

        def max_intron(x):
            array = x.split(",")
            if len(array) == 2:
                return 0
            else:
                intron_length = [
                    int(array[i]) - int(array[i - 1]) for i in range(1, len(array) - 1)
                ]
                return max(intron_length)

        """Parsing the pep fasta file and saving it in a df"""
        pep = SeqIO.parse(self.target_faa, "fasta")
        pep_df = pd.DataFrame(
            [[record.id, record] for record in pep], columns=["name", "seq"]
        )

        """Parsing """
        isoform_df = pd.read_csv(self.target_isoform, header=0, sep="\t")
        """parsing the bed file of target"""
        bed_df = pd.read_csv(
            self.target_bed,
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
            sep="\t",
        )
        bed_tmp = bed_df.copy()
        bed_tmp.cdsStarts.apply(lambda x: max_intron(x))
        bed_tmp["max_intron"] = bed_tmp.cdsStarts.apply(lambda x: max_intron(x))
        bed_all = pd.merge(bed_tmp, pep_df, how="inner", on="name")
        bed_all = pd.merge(
            isoform_df, bed_all, how="inner", left_on="TransID", right_on="name"
        )

        return bed_all

    def _miniprot_index(self):
        check_point = self.check / "miniprot_index.ok"
        if not check_point.exists():
            print(f"miniprot: Indexing the query file")
            cmd = f"{self.miniprot} -t {self.threads} -d {self.root / '.miniprot_index'} {self.query_f}"
            shell = run_cmd(cmd, self.work_dir)
            shell.on_parally()
            check_point.touch()
        else:
            print(
                f"The check_point have been found at {check_point.absolute()}, passing this step"
            )
        return self.root / ".miniprot_index"

    def _miniprot_alignment(self):
        """Saving the sequence of each gene in corresponding directory"""

        def save_fasta(x, directory):
            single_dir = directory / x.GeneID.iloc[0]
            single_dir.mkdir(exist_ok=True, parents=True)
            SeqIO.write(x["seq"], single_dir / ".pep.fa", "fasta")

        """Genering the miniprot mapping commands. In this study, the max intron size (parameter -G) were based on different gene"""

        def generate_cmd(x, directory, miniprot_path, index_path):
            gene_directory = directory / x.GeneID.iloc[0]
            max_intron = math.ceil(int(max(x.max_intron)) / 1000)
            if max_intron == 0:
                G_parameter = ""
            else:
                G_parameter = f"-G {max_intron}k"
            cmd = f"{miniprot_path} --gff-only -t 4 {G_parameter} {index_path} {gene_directory / '.pep.fa'} > {gene_directory / 'miniprot.gff'}"
            return cmd

        work_dir = self.work_dir
        check_point = self.check / "miniprot_mapping.ok"
        # Filtering the gene based on 'select' file, and making directory for each gene
        select = [i.strip() for i in open(self.select, "r")]
        target_select = self.target_df.copy()
        target_select = target_select[target_select.iloc[:, 0].isin(select)]

        if not check_point.exists():
            target_select.groupby("GeneID", as_index=False).apply(
                lambda x: save_fasta(x, work_dir)
            )
            cmds = (
                target_select.groupby("GeneID", as_index=False)
                .apply(
                    lambda x: generate_cmd(
                        x, work_dir, self.miniprot, self.miniprot_index
                    )
                )
                .iloc[:, 1]
            )
            shell = run_cmd(cmds, run_dir=self.work_dir)
            shell.on_parally(
                bar=True,
                desc="Miniprot",
                unit="gene",
                threads=self.threads,
                error="miniprot_mapping.log",
            )
            check_point.touch()
        else:
            print(
                f"The check point have been found at {check_point.absolute()}, passing this step"
            )

        return target_select

    def _transform_2bit_to_fa(self, twoBit, fa):
        twoBit = Path(twoBit)
        fa = Path(fa)
        if not fa.exists():
            cmd = f"{self.twoBitToFa} {twoBit} {fa}"
            shell = run_cmd(cmd, twoBit.parent)
            shell.on_parally()
        return fa

    def _obtain_2bit_info(self, twoBit):
        twoBit = Path(twoBit)
        cmd = f"{self.twoBitInfo} {twoBit} stdout"
        shell = run_cmd(cmd, run_dir=self.root)
        shell_out = shell.on_parally(out=True)["output"][0]

        # 不知为何最后会有一个None
        return {
            i[0]: i[1]
            for i in [j.split("\t") for j in shell_out.split("\n")]
            if len(i) == 2
        }

    def _extract_lastz_align_region(self):
        """Parsing the output miniprot to obtained the region of query for lastz"""

        def parse_miniprot_gff(geneid, directory, qplank, query_info):
            gene_dir = directory / geneid
            gff_df = pd.read_csv(
                gene_dir / "miniprot.gff",
                sep="\t",
                header=None,
                comment="#",
                names=[
                    "seqid",
                    "source",
                    "type",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "phase",
                    "attributes",
                ],
            )
            if len(gff_df) != 0:
                gff_df = gff_df[gff_df.type == "mRNA"]
                gff_df_range = gff_df.groupby("seqid", as_index=False).agg(
                    {"start": np.min, "end": np.max}
                )

                return gff_df_range.apply(
                    lambda x: [geneid] + overlap(x, qplank, query_info), axis=1
                ).tolist()
                # for _, data in gff_df_range.iterrows():
                # return [x] + overlap(data, qplank, query_info)
                #    return [x] + [np.nan, 0, 0]
            else:
                return [[geneid, np.nan, 0, 0]]

        def parse_target_gff(x, tplank, target_info):
            return [x["GeneID"]] + overlap(x, tplank, target_info)

        def overlap(x, plank, query_info):
            chr_start = 0
            chr_end = int(query_info[x["seqid"]])
            start = int(x["start"]) - int(plank)
            end = int(x["end"]) + int(plank)
            if start < chr_start:
                new_start = 1
            else:
                new_start = start
            if end > chr_end:
                new_end = chr_end
            else:
                new_end = end
            # return pd.Series([x.seqid, new_start, new_end], index = ['chr', 'start', 'end'])
            return [x.seqid, new_start, new_end]

        query_info = self._obtain_2bit_info(self.query)
        target_info = self._obtain_2bit_info(self.target)

        target = self.target_df.copy()

        process_bar = tqdm(
            total=len(set(target.GeneID.tolist())),
            desc=f"{'Update query':25}",
            unit="gene",
        )
        query_update_region = list()
        for geneid in set(target.GeneID.tolist()):
            query_update_region += parse_miniprot_gff(
                geneid, self.work_dir, self.qplank, query_info
            )
            process_bar.update(1)
        process_bar.close()

        """ Parsing the bed file of target to obtaine the region of target for lastz"""
        target_gff = pd.read_csv(
            self.target_gff,
            header=None,
            comment="#",
            names=[
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
            sep="\t",
            low_memory=False,
        )
        target_gff = target_gff[target_gff.type == "gene"]
        target_gff["GeneID"] = target_gff.attributes.apply(
            lambda x: re.match(r"ID=gene:[\w\d]+", x).group().replace("ID=gene:", "")
        )

        target_gff = target_gff[target_gff.GeneID.isin(target.GeneID.tolist())]
        target_update_region = []
        for _, data in target_gff.iterrows():
            target_update_region.append(
                parse_target_gff(data, self.tplank, target_info)
            )

        lastz_align_region = pd.merge(
            pd.DataFrame(
                target_update_region, columns=["GeneID", "chr", "start", "end"]
            ),
            pd.DataFrame(
                query_update_region, columns=["GeneID", "chr", "start", "end"]
            ),
            on="GeneID",
            how="inner",
            suffixes=["_target", "_query"],
        )

        return lastz_align_region

    def results(self):
        return self.lastz_align_region


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
        dest="work_dir",
        help="The directory was used to save the results",
    )
    args = parse.parse_args()
    return args


def main():
    args = make_args()
    step_one = extract_region(args)
    lastz_align_region = step_one.results()


if __name__ == "__main__":
    main()
