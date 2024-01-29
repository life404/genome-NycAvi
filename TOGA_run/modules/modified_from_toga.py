#! /usr/bin/python3

"""
This py file are modified from some module files of toga
"""

from pathlib import Path
import numpy as np
import pandas as pd
from collections import defaultdict

class overlap_select:
    def __init__(self, bed_buff, chain):
        self.bed_buff = bed_buff
        self.chain = chain
    
    def parse_bed(self, bed_buff):
        ranges = []
        for line in self.bed_buff.split('\n'):
            if line == '':
                continue
            ranges.extend(self.make_bed_ranges(line))
        return list(sorted(ranges, key=lambda x: x[0]))     
        
    def make_bed_ranges(self, bed_line):
        """Convert a bed line to a set of exon ranges."""
        line_info = bed_line.split("\t")
        
        # parse bed-12 line according to specification
        chrom_start = int(line_info[1])
        gene_name = line_info[3]
        
        # basically get exon absolute coordinates
        block_starts = [chrom_start + int(x) for x in line_info[11].split(",") if x != ""]
        block_sizes = [int(x) for x in line_info[10].split(",") if x != ""]
        blocks_num = int(line_info[9])
        block_ends = [block_starts[i] + block_sizes[i] for i in range(blocks_num)]
        
        #create and return (exon start, exon end, gene name) tuples
        genes = [gene_name for _ in range(blocks_num)]
        return list(zip(block_starts, block_ends, genes))
    
    def chain_reader(self):
        """Yield chain blocks one by one"""
        chain_data = self.chain.split('\n')
        chain_head = self.chain_data[0]
        del chain_data[0]
        
        # define start point
        progress = int(chain_head.split()[5])
        blocks_num = len(chain_data)
        
        for i in range(blocks_num):
            block_info = chain_data[i].split()
            if len(block_info) == 1:
                block_size = int(block_info[0])
                block_start = progress
                block_end = block_start + block_size
                yield block_start, block_end
                break
                
            block_size = int(block_info[0])
            dt = int(block_info[1])
            block_start = progress
            block_end = block_start + block_size
            progress = block_end + dt 
            yield block_start, block_end
    
    def intersect(chain_block, bed_block):
        return min(chain_block[1], bed_block[1]) - max(chain_block[0], bed_block[0])
    
    def overlap_select(self):
        ranges = self.parse_bed(self.bed)
        genes = [x[2] for x in ranges]
        
        bed_overlaps = {gene: 0 for gene in genes}
        bed_covered_times = defaultdict(set)
        
        chain_len = 0
        start_with = 0
        
        for block in self.chain_reader(self.chain):
            
            chain_len = block[1] - block[0]
            bed_num = 0
            FLAG = False
            FIRST = True
            bed_num = 0
            
            while True:
                if FIRST:
                    bed_num = start_with
                else:
                    bed_num += 1 
                
                if bed_num >= len(ranges):
                    break
                
                exon = ranges[bed_num]
                
                if block[1] < exon[0]:
                    # The chain end is smaller than first exon
                    break
                
                intersection = self.intersect(block, (exon[0], exon[1]))
                if intersection > 0:
                    if not FLAG:
                        start_with = bed_num
                        FLAG = True
                    
                    bed_overlaps[exon[2]] += intersection
                    bed_covered_times[exon[2]].add(exon[0])
                
                else:
                    if block[0] > block[1]:
                        continue
                    elif FLAG:
                        break
        
        return chain_len, bed_overlaps, bed_covered_times
                        
class chain_runner:
    def __init__(self, bed_file, geneids):
        self.bed = bed_file
        self.chain_id = 
        self.chain = self.extract_chain(chain_file, chain_dict, chain_id)
        
    def extract_chain(self, chain_file, chain_dict, chain_id):
        