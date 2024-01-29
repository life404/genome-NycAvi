#! /usr/bin/python3

from pathlib import Path
from Bio import SeqIO
import pandas as pd
import numpy as np
import shutil
import argparse
from tqdm import tqdm
import concurrent.futures as cf

from modules.faToTwoBit_cmd import faToTwoBit


def target_split(target, target_bed, run_dir, tplank, gene_lst, rmfa, threads=32):
    target = Path(target)
    target_bed = Path(target_bed)
    run_dir = Path(run_dir)
    Parts = run_dir / "Parts"
    Parts.mkdir(exist_ok=True, parents=True)

    target_tmp = target_split_by_chromosome(target=target)

    bed_df = pd.read_csv(
        target_bed,
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

    process_bar = tqdm(
        total=len(gene_lst), desc=f"{'Splitting target fasta':25}", unit="gene"
    )
    with cf.ProcessPoolExecutor(max_workers=threads) as e:
        process_lst = [
            e.submit(target_split_sub, gene, bed_df, tplank, target_tmp, Parts, rmfa)
            for gene in gene_lst
        ]
        for process in cf.as_completed(process_lst):
            process.result()
            process_bar.update(1)
    process_bar.close()

    return gene_lst


def target_split_sub(gene, bed_df, tplank, target_tmp, run_dir, rmfa):
    isoform = Path(run_dir) / f"{gene}" / "isoform.tsv"
    isoform_df = pd.read_csv(isoform, sep="\t", header=0)
    bed_df = bed_df[bed_df.name.isin(isoform_df.TransID)]
    region = bed_df.groupby("chrom", as_index=False).agg(
        {"mRNAStart": np.min, "mRNAEnd": np.max}
    )
    region_chr = region.chrom[0]
    region_start = int(region.mRNAStart[0])
    region_end = int(region.mRNAEnd[0])
    tplank_left = int(tplank)
    tplank_right = int(tplank)
    record = SeqIO.read(target_tmp / f"{region_chr}.fa", "fasta")

    chrStart = 0 if (region_start - tplank_left) <= 0 else region_start - tplank_left
    chrEnd = len(record.seq)
    chrEnd = (
        chrEnd if region_end + tplank_right >= chrEnd else region_end + tplank_right
    )

    target_fa = Path(run_dir) / f"{gene}" / "target.fa"
    SeqIO.write(record[chrStart:chrEnd], target_fa, "fasta")
    faToTwoBit(target_fa, target_fa.parent, rmfa=rmfa)

    bed_df.loc[:, "mRNAStart"] = bed_df.loc[:, "mRNAStart"] - chrStart
    bed_df.loc[:, "mRNAEnd"] = bed_df.loc[:, "mRNAEnd"] - chrStart
    bed_df.loc[:, "exonStart"] = bed_df.loc[:, "exonStart"] - chrStart
    bed_df.loc[:, "exonEnd"] = bed_df.loc[:, "exonEnd"] - chrStart
    modified_bed = Path(run_dir) / f"{gene}" / "toga.bed"
    bed_df.to_csv(modified_bed, sep="\t", index=False, header=False)


def target_split_rerun(dir1, dir2):
    dir1 = Path(dir1).absolute()
    dir2 = Path(dir2).absolute()
    path1 = dir1 / "target.2bit"
    path2 = dir2 / "target.2bit"
    bed1 = dir1 / "toga.bed"
    bed2 = dir2 / "toga.bed"
    isoform1 = dir1 / "isoform.tsv"
    isoform2 = dir2 / "isoform.tsv"

    shutil.copy(path1, path2)
    shutil.copy(bed1, bed2)
    shutil.copy(isoform1, isoform2)
    return path2

def target_split_by_chromosome(target):
    target = Path(target).absolute()
    target_split = target.parent / '.target.split'
    target_split.mkdir(exist_ok=True, parents=True)
    if list(target_split.glob('*.2bit')):
        pass
    else:
        target_records = SeqIO.to_dict(SeqIO.parse(target, "fasta"))
        for index, value in target_records.items():
            split_file = target_split / f'{index}.fa'
            SeqIO.write(value, split_file, "fasta")
            faToTwoBit(split_file, target_split)
        faToTwoBit(target, target_split)
    return target_split   


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-t", dest="target", help="target file")
    parse.add_argument("--bed", dest="bed", help="bed file")
    parse.add_argument("-o", dest="run_dir", help="working directory")
    parse.add_argument("-s", dest="select", help="selected gene")
    args = parse.parse_args()

    tplank = 5000
    gene_lst = [i.strip() for i in open(args.select, "r")]

    target_split(args.target, args.bed, args.run_dir, tplank, gene_lst)


if __name__ == "__main__":
    main()
