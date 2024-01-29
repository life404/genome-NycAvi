#! /usr/bin/python3

from pathlib import Path
import pandas as pd
import sys
from concurrent import futures as cf


def parse(run_dir, output, selected):
    run_dir = Path(run_dir).absolute()
    loss_info_list = list()
    for gene in selected:
        toga_info = run_dir / "Parts" / gene / "toga" / "loss_summ_data.tsv"
        if toga_info.exists():
            for line in open(toga_info, "r"):
                if line.startswith("GENE"):
                    loss_info_list.append(line.strip().split())
        else:
            continue
    loss_info_df = pd.DataFrame(
        loss_info_list, columns=["type", "gene", "status"], index=None
    )

    unknow_genes = list(set(selected).difference(set(loss_info_df.gene.tolist())))
    unknow_gene_df = pd.DataFrame(
        {
            "type": ["U"] * len(unknow_genes),
            "gene": unknow_genes,
            "status": ["E"] * len(unknow_genes),
        }
    )
    loss_info_df = pd.concat([loss_info_df, unknow_gene_df])
    loss_info_df.to_csv(Path(output).absolute(), header=True, index=False, sep="\t")


def main():
    run_dir = sys.argv[1]
    output = sys.argv[2]
    parse(run_dir, output)


if __name__ == "__main__":
    main()
