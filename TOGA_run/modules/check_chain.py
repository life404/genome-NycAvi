#! /usr/bin/python3

import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import sys
pd.set_option('display.max_rows', 1000)


def check_chain(dir_lst, score_threshold):
    blank_chain_lst = list()
    process_bar = tqdm(total=len(dir_lst), desc='checking', unit='gene')
    for directory in dir_lst:
        chain = Path(directory).absolute() / "lastz" / "clean.chain"
        file_buffer = [i.strip().split() for i in open(chain, 'r') if i.startswith('chain')]
        if len(file_buffer) == 0:
            blank_chain_lst.append(directory.stem)
        else:
            file_buffer = [i.strip().split() for i in open(chain, 'r') if i.startswith('chain')]
            df = pd.DataFrame(
                file_buffer,
                columns=[
                    "chain",
                    "score",
                    "tName",
                    "tSize",
                    "tStrand",
                    "tStart",
                    "tEnd",
                    "qName",
                    "qSize",
                    "qStrand",
                    "qStart",
                    "qEnd",
                    "id",
                ],
            )
            df.score = df.score.astype(int)
            score_max = int(np.max(df.score))
            if score_max < score_threshold:
                blank_chain_lst.append(directory.stem)
        process_bar.update(1)
    process_bar.close()
    return blank_chain_lst


def main():
    dir_lst = [Path(i.strip()).absolute() for i in open('/home/panda2bat/Avivorous_bat/output/15_gene_loss/DesRot/output2/checkpoints/lastz.parallel.temp', 'r')]
    score_threshold = 15000
    blank_chain_lst = check_chain(dir_lst, score_threshold)
    print(len(blank_chain_lst))
    
if __name__ == "__main__":
    main()
    
    

    
    
