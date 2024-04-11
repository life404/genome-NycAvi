#! /usr/bin/env python3

import argparse
from pathlib import Path
import subprocess


def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument("--genome", dest = "genome",type = str, default = '/home/panda2bat/Avivorous_bat/output/04_genome-polish/genome.nextpolish.fasta', help = "The path of polished fasta file")
    parse.add_argument("--hic_path", dest = "hicpath",type = str, default = '/home/panda2bat/Avivorous_bat/input/HiC/N.aviator', help = "The path of directory, which contains all hic fastq files")
    parse.add_argument("--ou_path", dest = "oupath", type = str, help = 'The path of output directory')
    parse.add_argument("--step", dest = "step", type = int, help = 'The step of programe, 1:chromap_index; 2:chromap_alignment; 3:yash_scaffold; 4:yash2juicer; 5:manual_curation_juicer')
    parse.add_argument("--review", dest='review', type=str, help='The review results of juicer, only be used in second round of 3D-DNA')
    args = parse.parse_args()
    return(args)

def cmd(command):
    run = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    run_stdout, run_stderr = run.communicate()
    if run.poll() == 0:
        for line in run_stdout:
            print(line.decode('utf-8', 'ignore'))
    else:
        for line in run_stderr:
            print('Something is wrong:\n')
            print(line)
            exit(1)

def chromap_index(genome, ou_path):
    genome_path = Path(ou_path) / 'reference.polised.fasta'
    
    if not genome_path.exists():
        genome_path.symlink_to(Path(genome))
    
    chromap_index_path = Path(ou_path) / 'reference.polised.index'
    
    samtools_index_command = f'/home/panda2bat/TOOLS/samtools-last/bin/samtools faidx {str(genome_path)}'
    chromap_index_command = f'/home/panda2bat/TOOLS/miniconda3/envs/hic-scaffolding/bin/chromap -i -r {str(genome_path)} -o {str(chromap_index_path)}'
    
    if not (Path(ou_path) / 'reference.polished.fastas.fai').exists():
        cmd(samtools_index_command)
        
    if not chromap_index_path.exists():
        cmd(chromap_index_command)
    
    return(genome_path, chromap_index_path)
    
def chromap_alignment(ou_path, hic_path, genome_path, chromap_index_path):
    sam_path = Path(ou_path) / 'chromap.aligned.sam'
    hic1_lst = [str(i) for i in list(Path(hic_path).glob('*.fq.gz')) if '_1.fq.gz' in i.name]
    hic2_lst = [str(i) for i in list(Path(hic_path).glob('*.fq.gz')) if '_2.fq.gz' in i.name]
    hic1_lst.sort()
    hic2_lst.sort()
    
    chromap_alignment_command = f"/home/panda2bat/TOOLS/miniconda3/envs/hic-scaffolding/bin/chromap --preset hic -r {str(genome_path)} -x {str(chromap_index_path)} --remove-pcr-duplicates -1 {','.join(hic1_lst)} -2 {','.join(hic2_lst)} --SAM -o {str(sam_path)} -t 250"
    sam2bam_command = f"/home/panda2bat/TOOLS/miniconda3/envs/hic-scaffolding/bin/samtools sort -@ 240 -O BAM -o {str(sam_path).replace('.sam', '.bam')} {str(sam_path)}"
    
    print(f"run command: {chromap_alignment_command}")
    if not sam_path.exists():
        cmd(chromap_alignment_command)
    print(f"finished")
    print("################################################")
    
    print(f"run command: {sam2bam_command}")
    if not Path(str(sam_path).replace('.sam', '.bam')).exists():
        cmd(sam2bam_command)
    print(f"finishde")

def yash_scaffold(ou_path):
    bam_path = Path(ou_path) / 'chromap.aligned.bam'
    genome_path = Path(ou_path) / 'reference.polised.fasta'
    output_path = Path(ou_path) / 'yash_scaffold.fasta'
    yash_command = f"/home/panda2bat/TOOLS/miniconda3/envs/hic-scaffolding/bin/yahs {genome_path} {bam_path} -o {output_path}"
    print(yash_command)
    cmd(yash_command)

def yash2juicer(ou_path):
    bin_path = Path(ou_path) / 'yash_scaffold.fasta.bin'
    gpg_path = Path(ou_path) / 'yash_scaffold.fasta_scaffolds_final.agp'
    fai_path = Path(ou_path) / 'reference.polised.fasta.fai'
    out_JBAT_path = Path(ou_path) / 'out_JBAT'
    yash2juicer_command = f'/home/panda2bat/TOOLS/miniconda3/envs/hic-scaffolding/bin/juicer pre -a -o {out_JBAT_path} {bin_path} {gpg_path} {fai_path} >{out_JBAT_path}.log 2>&1'
    print(yash2juicer_command)


def hic_D3(ou_path, genome):
    root_path = Path(ou_path) / "3d-dna"
    ou_path = root_path / 'first'
    ou_path.mkdir(exist_ok=True, parents=True)

    # The merged_nodups.txt should generate manually and be removed to the root directory
    
    D3_first_command = f"~/TOOLS/3d-dna/run-asm-pipeline.sh -r2 {genome} {root_path}/merged_nodups.txt"
    print(D3_first_command)

def hic_D3_second(ou_path, genome, review):
    root_path = Path(ou_path) / '3d-dna'
    ou_path = root_path / 'second'
    ou_path.mkdir(exist_ok=True, parents=True)

    D3_second_command = f"~/TOOLS/3d-dna/run-asm-pipeline.sh -r {review} {genome} {root_path}/merged_nodups.txt"
    print(D3_second_command)

def manual_curation_juicer(ou_path):
    out_JBAT_txt = Path(ou_path) / 'out_JBAT.txt'
    out_JBAT_hic = Path(ou_path) / 'out_JBAT.hic'
    out_JBAT_log = Path(ou_path) / 'out_JBAT.log'
    with open (out_JBAT_log, 'r') as file:
        for line in file:
            if line.startswith('PRE_C_SIZE'):
                genomeSize = line.split()[-1]
    
    juicer_tools_command = f'java -Xmx36G -jar /home/panda2bat/TOOLS/juicer-1.6/CPU/common/juicer_tools.1.9.9_jcuda.0.8.jar pre {out_JBAT_txt} {out_JBAT_hic} assembly {genomeSize}'
    print(juicer_tools_command)

def juicer_post(ou_path):
    review_assembly = Path(ou_path) / 'out_JBAT.review.assembly'
    liftover_agp = Path(ou_path) / 'out_JBAT.liftover.agp'
    genome = Path(ou_path) / 'reference.polised.fasta'
    juice_post_command = f'/home/panda2bat/TOOLS/miniconda3/envs/hic-scaffolding/bin/juicer post -o out_JBAT {review_assembly} {liftover_agp} {genome}'
    print(juice_post_command)

def hic_contact_map(ou_path):
    final_agp = Path(ou_path) / 'out_JBAT.FINAL.agp'
    yahs_bin = Path(ou_path) / 'yash_scaffold.fasta.bin'
    ref_index = Path(ou_path) / 'out_JBAT.FINAL.fa.fai'
    
    samtools_faidx = f"/home/panda2bat/TOOLS/miniconda3/envs/hic-scaffolding/bin/samtools faidx {str(ref_index).replace('.fa.fai', '.fa')}"
    if not ref_index.exists():
        cmd(samtools_faidx)
    
    juicer_pre = f"(/home/panda2bat/TOOLS/miniconda3/envs/hic-scaffolding/bin/juicer pre {yahs_bin} {final_agp} {ref_index} | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G|awk 'NF' > out_CONTACT.txt.part) && (mv out_CONTACT.txt.part out_CONTACT.txt)"
    juicertools_pre = f"java -Xmx36G -jar /home/panda2bat/TOOLS/juicer-1.6/CPU/common/juicer_tools.1.9.9_jcuda.0.8.jar pre out_CONTACT.txt out_CONTACT.hic <(cut -f1,2 {ref_index})"
    print(juicer_pre)
    print(juicertools_pre)
 
def main():
    args = make_parse()
    genome = args.genome
    hic_path = args.hicpath
    ou_path = args.oupath
    review = args.review
    step = args.step
    
    Path(ou_path).mkdir(exist_ok = True, parents = True)
    
    if step == 1:
        genome_path, chromap_index_path = chromap_index(genome, ou_path)
    if step == 2:
        chromap_alignment(ou_path, hic_path, genome_path, chromap_index_path)
    if step == 3:
        yash_scaffold(ou_path)
    if step == 4:
        yash2juicer(ou_path, genome)
    if step == 5:
        manual_curation_juicer(ou_path)
    if step == 6:
        juicer_post
    if step == 7:
        hic_contact_map(ou_path)
        
    
if __name__ == "__main__":
    main()
    
    
    
            
