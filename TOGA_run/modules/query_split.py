#! /usr/bin/python3

from pathlib import Path
from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
import concurrent.futures as cf
import time

from modules.miniprot_cmd import miniprot_index
from modules.miniprot_cmd import miniprot_mapping
from modules.parallel_execute import parallel_run_cmds, shell_cmd_run
from modules.faToTwoBit_cmd import faToTwoBit


def query_split(
    query,
    target_pep,
    isoform,
    gene_lst,
    run_dir,
    qplank=10 * 1000 * 1000,
    threads=32,
    error_path=None,
    check_path=None,
    rmfa=False,
):
    query = Path(query).absolute()
    target_pep = Path(target_pep).absolute()
    isoform = Path(isoform).absolute()
    run_dir = Path(run_dir).absolute()
    Parts = run_dir / "Parts"
    Parts.mkdir(exist_ok=True, parents=True)

    # preform miniprot index
    print("Indexing the query fasta")
    query_index = miniprot_index(query, run_dir, check_path=check_path)

    # paralle run miniprot mapping
    print("\n" + "Step1.5: miniprot mapping".center(70, "-") + "\n")
    mapping_check = run_dir / "checkpoints" / "miniprot_mapping.ok"

    if not mapping_check.exists():
        # merging the pep sequence with the isoform file using DataFrame
        pep_records = {}
        pep_records["id"] = [record.id for record in SeqIO.parse(target_pep, "fasta")]
        pep_records["seq"] = [
            str(record.seq) for record in SeqIO.parse(target_pep, "fasta")
        ]
        pep_dataframe = pd.DataFrame(pep_records)
        isoform_dataframe = pd.read_csv(isoform, sep="\t", header=0)
        isoform_dataframe = pd.merge(
            isoform_dataframe,
            pep_dataframe,
            how="inner",
            left_on="TransID",
            right_on="id",
        )

        # filter the isoform dataframe based on the selected gene list
        isoform_dataframe = isoform_dataframe[isoform_dataframe.GeneID.isin(gene_lst)]
        cmd_lst = isoform_dataframe.groupby("GeneID").apply(
            query_split_sub, run_dir=Parts, index_file=query_index
        )
        cmd_lst.to_csv(
            mapping_check.parent / "miniprot_mapping.temp",
            sep="\t",
            header=False,
            index=True,
        )

        process_bar = tqdm(
            total=len(cmd_lst.tolist()), desc=f"{'Miniprot mapping':25}", unit="gene"
        )
        with cf.ProcessPoolExecutor(max_workers=threads) as e:
            process_lst = [
                e.submit(shell_cmd_run, command=cmd, run_dir=run_dir)
                for cmd in cmd_lst.tolist()
            ]
            for process in cf.as_completed(process_lst):
                process.result()
                process_bar.update(1)
        process_bar.close()
        mapping_check.touch()

        # parallel_run_cmds(cmd_lst.tolist(), run_dir = run_dir, num = threads, error=f"minimap_mapping", error_path=error_path)
    else:
        print(f"The checkpoint file {mapping_check} had been found, pass this step")
        cmd_lst = pd.read_csv(
            mapping_check.parent / "miniprot_mapping.temp",
            sep="\t",
            header=None,
            index_col=0,
        )

    # extract region
    print("\n" + "Spliting the query fasta")
    query_tmp = query_split_by_chromosome(query)
    
    # Some genes may can't be found in query species
    final_lst = list()
    onlyTargetGene_lst = list()
    process_bar = tqdm(
        total=len(cmd_lst.index), desc=f"{'Spliting query fasta':25}", unit="gene"
    )
    with cf.ProcessPoolExecutor(max_workers=threads) as e:
        process_lst = [
            e.submit(
                query_mapped_region_extract,
                gene,
                Parts,
                qplank,
                query_tmp,
                error_path,
                rmfa,
            )
            for gene in cmd_lst.index
        ]
        for process in cf.as_completed(process_lst):
            flag, gene = process.result()
            if flag == "null":
                onlyTargetGene_lst.append(gene)
            final_lst.append(gene)
            process_bar.update(1)
    process_bar.close()

    # write these target specific genes in file
    if onlyTargetGene_lst:
        print(
            f"These are {len(onlyTargetGene_lst)} genes can't be identified in genome of query species.\nPlease check these genes in files {run_dir}/onlyTargetGene.txt."
        )
        onlyTargetGene_buff = open(run_dir / "onlyTargetGene.txt", "w")
        onlyTargetGene_buff.writelines([f"{i}\n" for i in onlyTargetGene_lst])
        onlyTargetGene_buff.close()

    return final_lst


# complex function of dataframe
def query_split_sub(x, run_dir, index_file):
    x = x.reset_index(drop=True)
    gene = x.GeneID[0]

    # generate working directory
    tmp_dir = Path(run_dir) / f"{gene}"
    tmp_dir.mkdir(exist_ok=True, parents=True)

    # write isoform file
    isoform = tmp_dir / "isoform.tsv"
    x.to_csv(
        isoform,
        sep="\t",
        header=True,
        index=False,
        columns=["GeneID", "TransID"],
    )

    # writing corresponding protein sequencs in workding directory
    pep_path = tmp_dir / ".pep.fa"
    pep_buff = open(pep_path, "w")
    for index, value in x.iterrows():
        pep_buff.write(f">{value.TransID}\n{value.seq}\n")
    pep_buff.close()

    # generating the miniprote mapping commands
    cmd = miniprot_mapping(pep_path, index_file, tmp_dir, onlycmd=True)
    return cmd


# First, extrac the mapped region based on gff file. Second, extend specific length of this region. Third, obtained the sequence, which will be used in lastz alignment
def query_mapped_region_extract(gene, Parts, qplank, query_tmp, error_path, rmfa):
    gff = Path(Parts) / f"{gene}" / "miniprot.gff"
    gff_df = pd.read_csv(
        gff,
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
    if not gff_df.empty:
        gff_df = gff_df[gff_df.type == "mRNA"]

        """
        if a gene have more than one transcript, the region should be the longest one
        Gene:             ----------------------------------
        Tran1:               ****************************
        Tran2:             ***************************
        Tran3:               ********************************
        region: -----------------------------------------------------------
                |   plank  |       mapping region             |   plank   |
        """
        # print(gff)
        gff_region = gff_df.groupby("seqid", as_index=False).agg(
            {"start": np.min, "end": np.max}
        )
        # if gene == 'ENSG00000115762':
        #    print(gff_region)
        # region_chr = gff_region.seqid[0]
        # region_start = int(gff_region.start[0])
        # region_end = int(gff_region.end[0])
        # region_start = region_start if region_start - int(qplank) > 0 else 0
        # region_end = int(region_end) + int(qplank)

        # query_fa = Path(Parts) / f"{gene}" / 'query.fa'
        # query_fa_buff = open(query_fa, 'w')
        # query_fa_buff.write(f">{region_chr}\n{str(query_records[region_chr][region_start:region_end].seq)}\n")
        # query_fa_buff.close()
        # return('ok')
        query_fa = Path(Parts) / f"{gene}" / "query.fa"
        record_lst = list()
        for index, value in gff_region.iterrows():
            region_chr = value.seqid
            region_start = int(value.start)
            region_end = int(value.end)
            region_start = (
                region_start - int(qplank) if region_start - int(qplank) > 0 else 0
            )
            region_end = region_end + int(qplank)

            record = SeqIO.read(query_tmp / f"{region_chr}.fa", "fasta")
            record_lst.append(record[region_start:region_end])
        SeqIO.write(record_lst, query_fa, "fasta")
        faToTwoBit(query_fa, query_fa.parent, rmfa=rmfa)
        return ("ok", gene)

    else:
        error_buff = open(error_path / f"miniprot_mapping_result.{gene}.error", "w")
        error_buff.write(
            f"The gene {gene} can be mapped to query genome.\nAs you can see, the gff file of miniprot is blank.\n{gene} will not be used in downstream analyssi.\nThis case may not a error. In sometimes, species may do not have specific genes. Such as, IL32 is only in some Simians.\nYou may need to check wheather this gene is in genome of query species"
        )
        error_buff.close()
        return ("null", gene)


def query_split_rerun(query, run_dir, rmfa=False):
    run_dir = Path(run_dir).absolute()
    query = Path(query).absolute()
    query_tmp = query.parent / ".query_split"
    gff = run_dir / "miniprot.gff"
    gff_df = pd.read_csv(gff, sep="\t", header=None, comment="#")
    chr_lst = set(gff_df.iloc[:, 0].to_list())
    query_fa = run_dir / "query.fa"
    SeqIO.write(
        [SeqIO.read(query_tmp / f"{chr}.fa", "fasta") for chr in chr_lst],
        query_fa,
        "fasta",
    )

    if rmfa:
        faToTwoBit(query_fa, run_dir, rmfa=rmfa)
    return run_dir / "query.fa"

def query_split_by_chromosome(query):
    query = Path(query).absolute()
    query_split = query.parent / ".query.split"
    query_split.mkdir(exist_ok=True, parents=True)
    if list(query_split.glob('*.2bit')):
        pass
    else:
        query_records = SeqIO.to_dict(SeqIO.parse(query, "fasta"))
        for index, value in query_records.items():
            split_file = query_split / f"{index}.fa"
            SeqIO.write(value, split_file, "fasta")
            faToTwoBit(split_file, query_split)
        faToTwoBit(query, query_split)
    return query_split
        
def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-q", dest="query", help="query file")
    parse.add_argument("-t", dest="target_pep", help="pep file of target")
    parse.add_argument("-i", dest="isoform", help="isoform file")
    parse.add_argument(
        "-g", dest="gene_lst", help="selected genes were used to in analysis"
    )
    parse.add_argument("-o", dest="run_dir", help="working directory")
    args = parse.parse_args()

    gene_lst = [i.strip() for i in open(args.gene_lst, "r")]

    query_split(
        query=args.query,
        target_pep=args.target_pep,
        isoform=args.isoform,
        gene_lst=gene_lst,
        run_dir=args.run_dir,
    )


if __name__ == "__main__":
    main()
