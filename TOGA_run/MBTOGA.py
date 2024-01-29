#! /usr/bin/python3

import argparse
from pathlib import Path

# sys.path.append(Path(__file__).resolve().parent)
from modules.extract_miniprot_alignment import extract_region
from modules.lastz_on_miniprot import lastz_align
from modules.toga_on_single_gene import run_toga


def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument(
        "-c", "--config", dest="config", type=str, help="The config file"
    )
    parse.add_argument(
        "-o",
        "--outdir",
        dest="workdir",
        type=str,
        help="The directory used to save the results",
    )
    parse.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=64,
        help="The number of threads",
    )
    parse.add_argument(
        "--stop_toga",
        dest='stop_toga',
        type=bool,
        default=False,
        help="Stop this programe before run toga"
    )
    args = parse.parse_args()
    return args


def main():
    args = make_parse()
    region_info = extract_region(args).run()
    region_info.to_csv(Path(args.workdir) / 'lastz.region', sep='\t', header=True, index=False)
    
    lastz = lastz_align(args)
    lastz_dir = lastz.run(region_info)
    
    if not args.stop_toga: 
        toga = run_toga(args)
        toga.results(lastz_dir)


if __name__ == "__main__":
    main()
