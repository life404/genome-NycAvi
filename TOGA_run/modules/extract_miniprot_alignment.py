#! /usr/bin/python3

from pathlib import Path
import numpy as np
import pandas as pd
import re
import argparse
from configparser import ConfigParser
from Bio import SeqIO
from modules.parallel_execute import run_cmd
from datetime import datetime

WORKDIR = "/home/panda2bat/Avivorous_bat/script/TOGA_run/test/output.new"
CONFIG = "/home/panda2bat/Avivorous_bat/script/TOGA_run/modules/config.ini"
THREADS = 128
SHIFT = 10
REVERSE = 10


class extract_region:
    def __init__(self, args) -> None:
        self.config = args.config
        
        # Path for results directories
        self.workdir = Path(args.workdir)
        self.workdir.mkdir(exist_ok=True, parents=True)
        self.checkpoint = self.workdir / 'checkpoint'
        self.checkpoint.mkdir(exist_ok=True, parents=True)
        self.log = self.workdir / 'log'
        self.log.mkdir(exist_ok=True, parents=True)
        
        # Path of software
        self.miniprot = self._parse_config('Software', 'miniprot')
        self.ucsc = self._parse_config('Software', 'ucsc_tools')
        self.twoBitToFa = self.ucsc / 'twoBitToFa'
        self.twoBitInfo = self.ucsc / 'twoBitInfo'
        
        # Path of input files
        self.query = self._parse_config('Path', 'query')
        self.target = self._parse_config('Path', 'target')
        self.bed = self._parse_config('Path', 'bed')
        # Read the isoform and pep file, the fasta in pep file but not in isoform file will be removed
        self.pep = self._parse_config('Path', 'faa')
        self.isoform = self._parse_config('Path', 'isoform')
        self.pep, self.isoform = self.filter_faa_on_isoform()
        
        # Other parameter
        self.threads = args.threads
        self.target_plank = self._parse_config('Parameter', 'tplank')
        self.query_plank = self._parse_config('Parameter', 'qplank')
        
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
    
    def filter_faa_on_isoform(self):
        # parsing isoform file，基于isoform筛选要研究的protein
        isoform_df = pd.read_csv(self.isoform, sep='\t', header=0)
        
        # filter pep file
        faa_filter_path = self.workdir / 'target.pep.fa'
        faa_df = pd.DataFrame([[record.id, record] for record in SeqIO.parse(self.pep, 'fasta')], columns=['TransID', 'record'])
        faa_df = faa_df[faa_df.TransID.isin(isoform_df.TransID.tolist())]
        SeqIO.write(faa_df.record.tolist(), faa_filter_path, 'fasta')
        
        return faa_filter_path, isoform_df
    
    def __twoBitToFa(self, twoBit):
        twoBit = Path(twoBit)
        fasta = twoBit.parent / f"{twoBit.stem}.fa"
        if fasta.exists():
            pass
        else:
            cmd = f"{self.twoBitToFa} {twoBit} {fasta}"
            execute = run_cmd(cmd, run_dir=twoBit.parent)
            execute.on_parally()
        return fasta
    
    def __twoBitInfo(self, twoBit):
        twoBit = Path(twoBit).resolve()
        info = twoBit.parent / f"{twoBit.stem}.info" 
        if not info.exists():
            cmd = f"{self.twoBitInfo} {twoBit} {info}"
            execute = run_cmd(cmd, run_dir=self.workdir)
            out = execute.on_parally()
        else:
            pass
        info_df = pd.read_csv(info, sep='\t', header=None, names=['chr', 'length'])
        return info_df
    
    def miniprot_run(self):
        index_path = self.__miniprot_index()
        gff_path = self.__miniprot_alignment(index_path)
        return gff_path
    
    def __miniprot_index(self):
        index_path = self.workdir / '.miniprot.index'
        query_fa = self.__twoBitToFa(self.query)
        cmd = f"{self.miniprot} -t {self.threads} -d {index_path} {query_fa}"
        
        if index_path.exists():
            pass
        else:
            print(f'Indexing the query genome ......')
            execute = run_cmd(cmd, run_dir=self.workdir)
            execute.on_parally()
        
        return index_path
    
    def __miniprot_alignment(self, index_path):
        checkpoint = self.workdir / 'checkpoint' / 'miniprot.ok'
        miniprot_dir = self.workdir / 'miniprot'
        miniprot_dir.mkdir(exist_ok=True, parents=True)
        gff_path = miniprot_dir / 'miniprot.gff'
        
        max_intron = int(self.__max_intron_size_on_bed(self.bed) / 1000)
        
        if checkpoint.exists():
            print('passing')
            pass
        else:
            cmd = f"{self.miniprot} --gff -t {self.threads} -G {max_intron}k {index_path} {self.pep} > {gff_path}"
            execute = run_cmd(cmd, run_dir=miniprot_dir)
            out = execute.on_parally()
            if sum(out['code']) == 0:
                checkpoint.touch()
            else:
                error_log = self.log / 'miniprot.error'
                error_buff = open(error_log, 'w')
                error_buff.writelines(out['output'])
                error_buff.close()
                print(f'The miniprot alignment have errors, please check the log {error_log}')
                exit(1)
        
        return gff_path
            
    def parse_gff_file(self, gff):
        gff_file = Path(gff)

        paf_df = pd.DataFrame(
            [i.strip().split() for i in open(gff_file, "r") if i.startswith("##PAF")]
        )
        gff_df = pd.read_csv(gff_file, header=None, sep="\t", comment="#")
        
        # analyzed the gff DataFrame
        gff_df = gff_df[[0, 2, 3, 4, 8]]
        gff_df.columns = ['chr', 'type', 'start', 'end', 'attributes']
        gff_df = gff_df[gff_df.type != 'stop_codon']
        gff_df['ID'] = gff_df.attributes.apply(lambda x: x.split(';')[0].split('=')[1])
        gff_df['gene'] = gff_df.attributes.apply(lambda x: re.search(r'Target=[\w\d]+', x).group().split('=')[1])
        exon_num = gff_df.groupby('ID', as_index=False).apply(lambda x: x.type.tolist().count('CDS'))
        exon_num.columns = ['ID', 'exon_num']
        gff_df = pd.merge(gff_df, exon_num, on='ID', how='left')
        gff_df = gff_df[gff_df.type == 'mRNA']
        gff_df = gff_df.drop(columns=['type', 'attributes'])
        

        # analyzed the paf df
        paf_df = paf_df[[1, 2, 6, 19]]
        paf_df.columns = ["gene", "pep_len", "chr", "cs"]
        paf_df["consensus"] = paf_df.cs.apply(
            lambda x: sum([int(i.strip(":")) for i in re.findall(r":\d+", x)])
        )
        paf_df["identity"] = paf_df.consensus.astype(int) / paf_df.pep_len.astype(int)
        paf_df = paf_df.drop(columns=['gene', 'cs', 'chr'])
        
        miniprot_info = pd.concat([gff_df.reset_index(drop=True), paf_df.reset_index(drop=True)], axis=1)
        
        return miniprot_info
    
    def parse_bed_file(self, bed):
        bed_df = pd.read_csv(bed, sep='\t', header=None)
        
        """ Sorting the gene position base on position """
        bed_df = bed_df.sort_values([0, 1])
        bed_df = bed_df[[0, 1, 2, 3, 9]]
        bed_df.columns = ['chr', 'start', 'end', 'gene', 'exon_num']
        return bed_df
    
    def __max_intron_size_on_bed(self, bed):
        bed_df = pd.read_csv(bed, sep='\t', header=None)
        
        """ Caculate the max size of intron """
        def intron_size(x):
            array = x.strip(',').split(',')
            if len(array) == 1:
                return 0
            else:
                return max([int(array[i]) - int(array[i-1]) for i in range(1, len(array))])
             
        return max(bed_df[11].apply(lambda x: intron_size(x)))
        
    
    def select_anchor_gene(self, miniprot_info, bed_info):
        merge_info = pd.merge(bed_info, miniprot_info, on='gene', how='inner', suffixes=['_bed', '_miniprot'])
        
        def exon_percent(x):
            if x['exon_num_bed'] == x['exon_num_miniprot'] == 1:
                percent = 1 / 1.5
            else:
                percent = x['exon_num_miniprot'] / x['exon_num_bed']    
            return 1 / percent if percent > 1 else percent
        
        merge_info['adjust.percent'] = merge_info.apply(lambda x: exon_percent(x) * x['percent'], axis=1)
    
    def integrated_miniprot_information(self, isoform_df, bed_df, miniprot_df):
       
        # 使用isoform筛选研究的基因，不在isoform范围内的基因剔除
        merge_info = pd.merge(isoform_df, bed_df, left_on='TransID', right_on='gene', how='left')
        merge_info = pd.merge(merge_info, miniprot_df, on='gene', how='left', suffixes=['_bed', '_miniprot'])
        merge_info = merge_info.sort_values(['chr_bed', 'start_bed']).reset_index(drop=True)
        
        return merge_info
        

    def adjust_miniprot_results(self, merge_info):

        # using the exon numeric to normalized the identity
        def normaliz_factor(x):
            # We punished the single exon transcripts.
            if x['exon_num_bed'] == ['exon_num_miniprot'] == 1:
                factor =  1 / 1.5
            else:
                factor = x['exon_num_miniprot'] / x['exon_num_bed']
            return (1 / factor) if factor > 1 else factor
        
        def adjust_and_update(x, passed):
            chr_bed = x['chr_bed']
            chr_miniprot = x['chr_miniprot']
            start = int(x['start_bed'])
            end = int(x['end_bed'])
            shift = passed[(passed.start_bed < start) & (passed.chr_bed == chr_bed)].index.tolist()
            reverse = passed[(passed.start_bed > start) & (passed.chr_bed == chr_bed)].index.tolist()
            update_chr = set(passed.loc[shift[- SHIFT:], 'chr_miniprot'].tolist() + passed.loc[reverse[:REVERSE], 'chr_miniprot'].tolist())
            
            update_chr = [i for i in update_chr if i != chr_miniprot]
            
            # 对于比对长度过短的序列，miniprot 确定的位置可能存在错误，因此，除了对miniprot 找出的候选位置位置比对外，还需要对整条染色体进行比对。该长度阈值暂定为 20氨基酸
            if x['consensus'] <= 20:
                update_chr.append(chr_miniprot)
            
            return ','.join(update_chr) if update_chr else pd.NA
        
        def adjust_and_update_na(x, passed):
            chr_bed = x['chr_bed']
            chr1 = passed[passed.chr_bed == chr_bed].chr_miniprot.tolist()
            chr2 = list()
            for i in passed[passed.chr_bed == chr_bed].chr_miniprot.tolist():
                chr2.extend(i.split(','))
            update_chr = list(set(chr1 + chr2))
            return ','.join(update_chr)
                
        merge_info_na = merge_info[merge_info.chr_miniprot.isna()]
        merge_info = merge_info[~merge_info.chr_miniprot.isna()].copy()            
        merge_info.loc[:, 'adjust_identity'] = merge_info.apply(lambda x: normaliz_factor(x) * x['identity'], axis = 1)
        
        # The 0.75 were be used the threshold. If a transcripts with a adjust_identity less than the threshold, it will be check. 
        passed = merge_info[merge_info.adjust_identity > 0.75].reset_index(drop=True)
        tester = merge_info[merge_info.adjust_identity <= 0.75].reset_index(drop=True).copy()
        tester.loc[:, 'update_chr'] = tester.apply(lambda x: adjust_and_update(x, passed), axis = 1)
        
        merge_update = pd.concat([passed, tester], axis=0, join='outer').sort_values(['chr_bed', 'start_bed'])
        
        # If a gene can identity any region
        tester = merge_info_na.copy()
        tester.loc[:, 'update_chr'] = tester.apply(lambda x: adjust_and_update_na(x, merge_update), axis=1)
        merge_update = pd.concat([merge_update, tester], axis=0, join='outer').sort_values(['chr_bed', 'start_bed'])
        
        return merge_update
    
    def extrac_lastz_region(self, merge_info):
        
        def extrac_region(x, query_2bit_df, target_2bit_df, query_plank, target_plank):
            df = x
            target_region = df.groupby('chr_bed', as_index=False).agg({'start_bed': np.min, 'end_bed': np.max})
    
            query_region = df.groupby('chr_miniprot', as_index=False).agg({'start_miniprot': np.min, 'end_miniprot': np.max})
            target_update = target_region.apply(lambda x: update_region(x, target_2bit_df, target_plank), axis = 1)
            query_update = query_region.apply(lambda x: update_region(x, query_2bit_df, query_plank), axis = 1)
            
            update_chr = df[~df.update_chr.isna()]['update_chr'].tolist()
            if update_chr:
                chr_lst = list()
                for i in update_chr:
                    chr_lst.extend(i.split(','))
                chr_lst = list(set(chr_lst))
                update_chr_df = pd.DataFrame({'chr_miniprot': chr_lst, 'start_miniprot':[np.nan]*len(chr_lst), 'end_miniprot':[np.nan]*len(chr_lst)})
                query_update = pd.concat([query_update, update_chr_df])
            
            target_update_transform = target_update.loc[target_update.index.repeat(len(query_update))]
            target_update_transform = target_update_transform.reset_index(drop=True)
            query_update_transform = pd.concat([query_update]*len(target_update))
            query_update_transform = query_update_transform.reset_index(drop=True)
            
            return pd.concat([target_update_transform, query_update_transform], axis=1)
            
        def update_region(x, twoBitInfo, plank):
            start = x.iloc[1]
            end = x.iloc[2]
            chrom = x.iloc[0]
            length = twoBitInfo[twoBitInfo.chr == chrom]['length'].iloc[0]
            new_s = 1 if int(start) - int(plank) <= 0 else int(start) - int(plank)
            new_e = length if int(end) + int(plank) > length else int(end) + int(plank)
            x.iloc[1] = new_s
            x.iloc[2] = new_e
            return x
            
        
        query_2bit_df = self.__twoBitInfo(self.query)
        target_2bit_df = self.__twoBitInfo(self.target)
        
        lastz_align_region = merge_info.groupby('GeneID').apply(lambda x: extrac_region(x, query_2bit_df, target_2bit_df, self.query_plank, self.target_plank)).reset_index()
        return lastz_align_region
    
    def run(self):

        print('Runing the miniprot alignment ......')
        start_time = datetime.now()
        gff_path = self.miniprot_run()
        end_time = datetime.now()
        print(f"{'time used':20}:{end_time - start_time}")

        print('Parsing the output of miniprot ......')
        start_time = datetime.now()
        miniprot_info = self.parse_gff_file(gff_path)
        bed_info = self.parse_bed_file(self.bed)
        merge_info = self.integrated_miniprot_information(self.isoform, bed_info, miniprot_info)
        print('Obtaining the region used for lastz ......')
        
        merge_update = self.adjust_miniprot_results(merge_info)
        lastz_align_region = self.extrac_lastz_region(merge_update)
        end_time = datetime.now()
        print(f"{'time used':20}:{end_time - start_time}")
        return lastz_align_region


def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument("--output", dest='workdir', type=str, default=WORKDIR, help='work directory')
    parse.add_argument("--config", dest="config", type=str, default=CONFIG, help="config file")
    parse.add_argument("--threads", dest="threads", type=str, default=THREADS, help="threads num")
    args = parse.parse_args()
    return args


def main():
    args = make_parse()
    lastz_align_region = extract_region(args).run()
    lastz_align_region.to_csv('/home/panda2bat/Avivorous_bat/script/TOGA_run/test/output.new/lastz_region', sep='\t', header=True, index=False)
    print(lastz_align_region)


if __name__ == "__main__":
    main()
