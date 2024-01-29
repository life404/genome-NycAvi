#! /usr/bin/python3

from pathlib import Path
import concurrent.futures as cf
from tqdm import tqdm
import argparse
from datetime import datetime
import psutil
import subprocess

from modules.query_split import (
    query_split,
    query_split_rerun,
    query_split_by_chromosome,
)
from modules.target_split import target_split, target_split_rerun
from modules.lastz_pipeline import lastz_pipeline, lastz_align_direct
from modules.lastz_rerun import lastz_rerun_on_cpu_percent
from modules.check_chain import check_chain
from modules.toga_cmd import toga, toga_base_on_memory
from modules.parse_toga import parse
from modules.faToTwoBit_cmd import faToTwoBit


def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument("-t", "--target", dest="target", help="The target fasta file")
    parse.add_argument("-q", "--query", dest="query", help="The query fasta file")
    parse.add_argument(
        "--bed", dest="bed", help="The bed12 format bed annotaion of target"
    )
    parse.add_argument(
        "-p", "--pep", dest="pep", help="The protein sequences fasta file of target"
    )
    parse.add_argument(
        "--isoform", dest="isoform", help="The isoform file of target file"
    )
    parse.add_argument(
        "-s", "--selected", dest="selected", help="The selected genes will be analyzed"
    )
    parse.add_argument("-o", "--output", dest="run_dir", help="The output directory")
    parse.add_argument(
        "--tplank", dest="tplank", type=int, default=20 * 1000, help="The tplank of target"
    )
    parse.add_argument(
        "--qplank",
        dest="qplank",
        type=int,
        default=10 * 1000 * 1000,
        help="The plank of query",
    )
    parse.add_argument(
        "--score_threshold",
        dest="score_threshold",
        type=float,
        default=15000,
        help="The identity threshold of miniprot mapping, the range is [0, 1]",
    )
    parse.add_argument(
        "--threads",
        type=int,
        default=32,
        dest="threads",
        help="The number of thread. In this programe the threads parameter set the number of thread used in miniprot and lastz. But the threads of toga will base on the aviliable memmory of machine and run memmory",
    )
    parse.add_argument(
        "--mem_max",
        type=int,
        dest="mem_max",
        default=300,
        help="The max memroy used in cesear process in toga.",
    )
    parse.add_argument(
        "--mem_min", type=int, dest="mem_min", default=5, help="memory min"
    )
    parse.add_argument(
        "--mem_step",
        type=int,
        dest="mem_step",
        default=20,
        help="The increased step value of memmory. If a gene need more memmory, the gene will be run with more memmory (previous + step)",
    )
    parse.add_argument(
        "--rm",
        default=True,
        type=bool,
        dest="rmfa",
        help="The intermidate fasta files will be removed to save the space",
    )

    args = parse.parse_args()
    return args


def main():
    args = make_parse()

    error_path = Path(args.run_dir) / "errors"
    error_path.mkdir(exist_ok=True, parents=True)

    check_path = Path(args.run_dir) / "checkpoints"
    check_path.mkdir(exist_ok=True, parents=True)

    selected_gene = [i.strip() for i in open(args.selected, "r")]
    
    # check the ulimit, if the ulimit set with the default value 1024, the downstram analysis may pull a error ' OSError: [Errno 24] Too many open files'
    ulimit_num = subprocess.Popen('ulimit -n', shell=True, stdout=subprocess.PIPE, stderr = subprocess.STDOUT, text=True).communicate()
    ulimit_num = int(ulimit_num[0].strip())
    if ulimit_num == 1024:
        print(f"The ulimit num is {ulimit_num}, please set with a bigger value")
        exit(1)

    # In this pipeline, the miniprot will be used. First, the protein sequences from target species will be aligned to genome file of query species. Then, the sequences in mapped region will be extracted from genome of query species. We also extrend the mapped region with a plank
    check = check_path / "query.split.ok"
    print("Step1: split query".center(70, "-") + "\n")
    if not check.exists():
        start_time = datetime.now()
        gene_lst = query_split(
            query=args.query,
            target_pep=args.pep,
            isoform=args.isoform,
            gene_lst=selected_gene,
            run_dir=args.run_dir,
            qplank=args.qplank,
            error_path=error_path,
            check_path=check_path,
            threads=args.threads,
            rmfa=args.rmfa,
        )
        check.touch()

        temp = check_path / "query.split.temp"
        temp_buff = open(temp, "w")
        temp_buff.writelines(f"{i}\n" for i in gene_lst)
        temp_buff.close()
        end_time = datetime.now()
        print(
            f"{'Start':<5}: {start_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'End':<5}: {end_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'Time':<5}: {end_time - start_time}"
        )
    else:
        print(f"The check point file {check} has beed found, pass this step")
        temp = check_path / "query.split.temp"
        gene_lst = [i.strip() for i in open(temp, "r")]

    # Be simililar with query, the coding region of each gene also be extracted from genom of target specie based on the bed annotation.
    check = check_path / "target.split.ok"
    print("\n" + "Step2: split target".center(70, "-") + "\n")
    if not check.exists():
        start_time = datetime.now()
        gene_lst = target_split(
            target=args.target,
            target_bed=args.bed,
            run_dir=args.run_dir,
            tplank=args.tplank,
            gene_lst=gene_lst,
            rmfa=args.rmfa,
            threads=args.threads,
        )
        check = check_path / "target.split.ok"
        check.touch()

        temp = check_path / "target.split.temp"
        temp_buff = open(temp, "w")
        temp_buff.writelines([f"{i}\n" for i in gene_lst])
        temp_buff.close()
        end_time = datetime.now()
        print(
            f"{'Start':<5}: {start_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'End':<5}: {end_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'Time':<5}: {end_time - start_time}"
        )
    else:
        print(f"The check point file {check} has been found, pass this step")
        temp = check_path / "target.split.temp"
        gene_lst = [i.strip() for i in open(temp, "r")]

    # run lastz in parallel
    # A gene without miniprot results will be passed in first round lastz pipeline
    miniprot_blank_lst = [
        i.strip() for i in open(Path(args.run_dir) / "onlyTargetGene.txt", "r")
    ]
    gene_lst = list(set(gene_lst).difference(set(miniprot_blank_lst)))
    check = check_path / "lastz.parallel.ok"
    print("\n" + "Step3: run lastz in parallel".center(70, "-") + "\n")
    if not check.exists():
        start_time = datetime.now()
        process_bar = tqdm(total=len(gene_lst), desc=f"{'lastz':25}", unit="gene")
        with cf.ProcessPoolExecutor(max_workers=args.threads) as e:
            process_lst = [
                e.submit(lastz_pipeline, gene, Path(args.run_dir) / "Parts", error_path)
                for gene in gene_lst
            ]
            lastz_out_lst = []
            for process in cf.as_completed(process_lst):
                out = process.result()
                lastz_out_lst.append(out)
                process_bar.update(1)
        process_bar.close()
        check.touch()
        temp = check_path / "lastz.parallel.temp"
        temp_buff = open(temp, "w")
        temp_buff.writelines([f"{i}\n" for i in lastz_out_lst])
        temp_buff.close()
        end_time = datetime.now()
        print(
            f"{'Start':<5}: {start_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'End':<5}: {end_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'Time':<5}: {end_time - start_time}"
        )
    else:
        print(f"The check point file {check} has been found, pass this step")
        temp = check_path / "lastz.parallel.temp"
        lastz_out_lst = [Path(i.strip()) for i in open(temp, "r")]

    # if rerun is true, the query split will re-extract suqences based on chrosome, not mapped region. Then the re-extract query sequences will be aligned by lastz
    check = check_path / "lastz.rerun.ok"
    print(
        "\n" + "Step 4: check and re-align the results of lastz".center(70, "-") + "\n"
    )
    if not check.exists():
        print(
            f"In this step, these genes will be need to re-align using lastz\n1. Miniprot can not align.\n2. The chain of lastz is blank\n3. The max score of chain blocks less than {int(args.score_threshold)}"
        )
        start_time = datetime.now()
        # The gene with blank chain, low quality (the max score of chain blocks less than score threshold) chains, or blank miniprot result will be aligned against to whole genome of query species by using lastz
        blank_chain_lst = check_chain(lastz_out_lst, int(args.score_threshold))
        blank_chain_lst += miniprot_blank_lst
        blank_chain_lst = list(set(blank_chain_lst))
        
        blank_chain_buff = open(Path(args.run_dir).absolute() / 'rerun_genes_lst.txt', 'w')
        blank_chain_buff.writelines([f"{i}\n" for i in blank_chain_lst])
        blank_chain_buff.close()
        #blank_chain_lst = [
        #    i.strip() for i in open(Path(args.run_dir) / "rerun_genes_lst.txt", "r")
        #]
        print(
            f"A total of {len(blank_chain_lst)} genes need to be re-aligned, these genes can be check on {Path(args.run_dir).absolute() / 'rerun_genes_lst.txt'}"
        )

        lastz_align_again_lst = lastz_rerun_on_cpu_percent(
            gene_lst=blank_chain_lst,
            query=args.query,
            root_dir=((Path(args.run_dir) / "Parts")),
            error_path=error_path,
        )

        lastz_out_lst = list(set(lastz_out_lst + lastz_align_again_lst))
        temp = check_path / "lastz_rerun.temp"
        temp_buff = open(temp, "w")
        temp_buff.writelines([f"{i}\n" for i in lastz_out_lst])
        temp_buff.close()
        check.touch()
        end_time = datetime.now()
        print(
            f"{'Start':<5}: {start_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'End':<5}: {end_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'Time':<5}: {end_time - start_time}"
        )
    else:
        print(f"The check point file {check} has been found, pass this step")
        temp = check_path / "lastz_rerun.temp"
        lastz_out_lst = [Path(i.strip()) for i in open(temp, "r")]

    # run toga in parallel
    lastz_out_lst = [Path(i.strip()) for i in open(temp, "r")]
    check = check_path / "toga.ok"
    print("\n" + "Step5: Run TOGA in parallel".center(70, "-") + "\n")
    if not check.exists():
        start_time = datetime.now()
        lastz_out_lst = sorted(lastz_out_lst, key=lambda x: x.stat().st_size)
        toga_run_dir = lastz_out_lst
        toga_finish_dir = []
        mem_lst = [args.mem_min] + list(
            range(0, (args.mem_max + args.mem_step), args.mem_step)[1:]
        )
        for cmx in mem_lst:
            if toga_run_dir:
                sub_start_time = datetime.now()
                if cmx == args.mem_min:
                    print("Run toga in first time")
                    cjn = 4
                else:
                    print(f"Rerun toga in {cmx} G memory")
                    cjn = 2
                mem = psutil.virtual_memory()
                available_mem = mem.available / 1024**3
                print(
                    f"{'Avialabel memory (G):':>26} {int(available_mem):<5}\n{'Memory limit of cesar (G):':>26} {cmx:<5}"
                )
                toga_finish_lst, toga_rerun_lst, toga_error_lst = toga_base_on_memory(
                    run_dir_lst=toga_run_dir, error_path=error_path, cmx=cmx, cjn=cjn
                )

                print(
                    f"A total of {len(toga_finish_lst)} genes have finished toga analysi"
                )
                print(
                    f"A total of {len(toga_rerun_lst)} genes need to run in more memory"
                )
                print(
                    f"A total of {len(toga_error_lst)} genes have error, the detailed error information can be checked in {error_path}/toga.*.error"
                )

                toga_run_dir = toga_rerun_lst
                toga_finish_dir = toga_finish_dir + toga_finish_lst
                sub_end_time = datetime.now()
                print(
                f"{'Time':<5}: {sub_end_time - sub_start_time}"
                )

        toga_temp_buff = open(check.parent / "toga.temp", "w")
        toga_temp_buff.writelines([f"{i}\n" for i in toga_finish_dir])
        toga_temp_buff.close()
        #toga_results_lst = toga_finish_dir

        check.touch()
        end_time = datetime.now()
        print(
            f"{'Start':<5}: {start_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'End':<5}: {end_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'Time':<5}: {end_time - start_time}"
        )
    else:
        print(f"The check point file {check} has been found, pass this step")
        #toga_results_lst = [i.strip() for i in open(check.parent / "toga.temp", "r")]

    print("\n" + "Step6: Parse the results of TOGA".center(70, "-") + "\n")
    start_time = datetime.now()
    parse_out = Path(args.run_dir).absolute() / "toga_status.info"
    parse(run_dir=args.run_dir, output=parse_out, selected=selected_gene)
    end_time = datetime.now()
    print(
        f"{'Start':<5}: {start_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'End':<5}: {end_time.strftime('%Y.%m.%d %H:%M:%S')}\n{'Time':<5}: {end_time - start_time}"
    )


if __name__ == "__main__":
    main()
