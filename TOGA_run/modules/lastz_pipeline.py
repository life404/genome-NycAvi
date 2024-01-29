#! /usr/bin/python3

from pathlib import Path
import concurrent.futures as cf
import sys
import shutil

from modules.lastz_cmd import lastz
from modules.lastz_cmd import axtChain
from modules.lastz_cmd import RepeatFiller
from modules.chainCleaner_cmd import chainCleaner
from modules.chainMergeSort_cmd import chainMergeSort


def lastz_pipeline(gene, run_dir, error=None, error_path=None):
    if gene:
        run_dir = Path(run_dir).absolute() / f"{gene}"
    else:
        run_dir = Path(run_dir).absolute()
        gene = run_dir.stem
    if not error:
        error = gene
    query_2bit = run_dir / "query.2bit"
    target_2bit = run_dir / "target.2bit"
    lastz_out = run_dir / "lastz"
    lastz_out.mkdir(exist_ok=True, parents=True)

    # preform lastz
    axt = lastz(
        target=target_2bit,
        query=query_2bit,
        run_dir=lastz_out,
        error=f"lastz.{error}",
        error_path=error_path,
    )

    # preform axtChain
    chain = axtChain(
        axt=axt,
        target=target_2bit,
        query=query_2bit,
        run_dir=lastz_out,
        error=f"axtChain.{error}",
        error_path=error_path,
    )

    # preform RepeatFiller
    filled_chain = RepeatFiller(
        chain=chain,
        target=target_2bit,
        query=query_2bit,
        run_dir=lastz_out,
        error=f"RepeatFiller.{error}",
        error_path=error_path,
    )
    
    
    # preform chainCleaner
    clean_chain = chainCleaner(
        chain=filled_chain,
        target=target_2bit,
        query=query_2bit,
        run_dir=lastz_out,
        error=f"chainCleaner.{error}",
        error_path=error_path,
    )

    return run_dir


def lastz_align_direct(
    gene, query_2bit, run_dir, error=None, error_path=None, threads=32
):
    run_dir = Path(run_dir).absolute() / f"{gene}"
    target_2bit = run_dir / "target.2bit"
    query_2bit = Path(query_2bit).absolute()
    query_2bit_slink = run_dir / "query.2bit"
    query_2bit_slink.symlink_to(query_2bit)
    lastz_out = run_dir / "lastz"
    lastz_out.mkdir(exist_ok=True, parents=True)

#    axt = lastz(
#        target=target_2bit,
#        query=query_2bit,
#        run_dir=lastz_out,
#        error=f"lastz.align_direct.{gene}",
#        error_path=error_path,
#    )
#    chain = axtChain(
#        axt=axt,
#        target=target_2bit,
#        query=query_2bit,
#        run_dir=lastz_out,
#        error=f"axtChain.align_direct.{gene}",
#        error_path=error_path,
#    )
#
#    filled_chain = RepeatFiller(
#        chain=chain,
#        target=target_2bit,
#        query=query_2bit,
#        run_dir=lastz_out,
#        error=f"RepeatFiller.align_direct.{gene}",
#        error_path=error_path,
#    )
#
#    if filled_chain.stat().st_size != 0:
#        clean_chain = chainCleaner(
#            chain=filled_chain,
#            target=target_2bit,
#            query=query_2bit,
#            run_dir=lastz_out,
#            error=f"chainCleaner.align_direct.{gene}",
#            error_path=error_path,
#        )
#        if clean_chain.stat().st_size == 0:
#            clean_chain = shutil.copy(filled_chain, lastz_out / "clean.chain")
#    else:
#        clean_chain = shutil.copy(chain, lastz_out / "clean.chain")
#
    return run_dir


def main():
    query = sys.argv[1]
    run_dir = sys.argv[2]
    query_2bit_split = list(
        Path(
            "/home/panda2bat/Avivorous_bat/output/15_gene_loss/DesRot/input/.query_split"
        ).glob("*.2bit")
    )

    lastz_align_direct(
        gene="ENSG00000126549",
        query=query,
        query_2bit_split=query_2bit_split,
        run_dir=run_dir,
    )


if __name__ == "__main__":
    main()
