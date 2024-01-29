#! /usr/bin/env python3

import argparse
from pathlib import Path
from utils import *
import random
import string
import re
import shutil
import pandas as pd
import numpy as np


def make_parse():
    parse = argparse.ArgumentParser()
    parse.add_argument(
        "-g", "--genome", dest="genome", help="The path of genome input (softmasked)"
    )
    parse.add_argument(
        "-p",
        "--protein",
        dest="protein",
        help="The path of protein fasta file of close relative, this protein file can be used in miniprot or braker3",
    )
    parse.add_argument(
        "-o", "--ou_path", dest="ou_path", help="The path of output directory"
    )
    parse.add_argument("-n", "--name", dest="prefix", help="The prefix of miniprot output file")
    parse.add_argument(
        "--rna",
        dest="rna",
        type=str,
        help="The path of directory, which contains RNA fastq files",
    )
    parse.add_argument(
        "--SS",
        dest="strand",
        type=str,
        default="false",
        choices=["false", "RF", "FR"],
        help="The BOOL value of SS_lib_type parameter in Trintiy",
    )
    parse.add_argument(
        "--dbname",
        dest = "dbname",
        type = str,
        help = 'The SQL name used in pasa, if not gived, the pasa db will named with 5 random letters',
        required=False
    )
    parse.add_argument(
        "--gff",
        dest = "gff",
        type = str, 
        help = "The GFF annotation of PASA update, this option only used in step 8 functional annotation."
    )
    parse.add_argument(
        "-s",
        "--step",
        dest="step",
        type=int,
        help="The step of annotation. 1:miniprot; 2:hisat2; 3:stringtie; 4:trinity+pasa; 5:braker; 6:EvidenceModeler; 7:PASA upate; 8:functional annotation"
    )
    parse.add_argument(
        "--verbose",
        dest="verbose",
        default=False,
        type=bool,
        help="Print the verbpse information (True) or not (Flase), default = False",
    )
    parse.add_argument(
        "--threads",
        dest = "threads",
        default = 32,
        type = int,
        help = "The number of threads"
    )
    args = parse.parse_args()
    return args

def miniprot_run(genome, ou_path, protein, prefix):
    ou_path = Path(ou_path) / "miniprot"
    ou_path.mkdir(exist_ok=True, parents=True)

    index_path = ou_path / f"{prefix}_index.mpi"
    paf_path = ou_path / f"{prefix}.paf"
    gff_path = ou_path / f"{prefix}.gff"
    aln_path = ou_path / f"{prefix}.aln"
    un_path = ou_path / f"{prefix}.unmapped.query.fa"
    gtf_path = ou_path / f"{prefix}.gtf"

    miniprot_index_command = f"/home/panda2bat/TOOLS/miniprot/miniprot-0.11_x64-linux/miniprot -t 240 -d {index_path} {genome}"
    miniprot_align_command = f"/home/panda2bat/TOOLS/miniprot/miniprot-0.11_x64-linux/miniprot -t 240 -I -u --gff {index_path} {protein} > {gff_path}"

    check = cmd_run_multiple(miniprot_index_command, ou_path, check="miniport_index")
    if check.exists():
        check = cmd_run_multiple(
            [miniprot_align_command], ou_path, check="miniport_align"
        )
    else:
        print(f"CAN NOT FOUND CHECK `{check}` OF PREVIOUS STEP, SOMETHING HAS WRONG")
        exit(1)
        
def sam2bam(ou_path):
    sam_lst = list(ou_path.glob("*.sam"))
    command_lst = []
    for sam in sam_lst:
        bam_path = ou_path / sam.name.replace(".sam", ".sort.bam")
        sam2bam_command = f"~/TOOLS/samtools-last/bin/samtools sort -@ 32 -O BAM -o {bam_path} {sam}\n"
        command_lst.append(sam2bam_command)

    check = cmd_run_multiple(command_lst, ou_path, check="sam2bam", num=8)
    return check

def bam_merge(ou_path):
    bam_lst = " ".join([str(i) for i in list(Path(ou_path).glob("*.bam"))])
    merge_path = Path(ou_path) / "merge.bam"
    samtools_merge = (
        f"~/TOOLS/samtools-last/bin/samtools merge -@ 128 {merge_path} {bam_lst}"
    )
    check = cmd_run_multiple([samtools_merge], ou_path, check="samtools_merge")

    return (check, merge_path)

def RNA_anno(rna_path, ou_path, genome, verbose):
    ou_path = Path(ou_path) / "RNA" / "hisat"
    ou_path.mkdir(exist_ok=True, parents=True)

    # build the index of hisat2
    hisat2_index = Path(ou_path) / "hisat2.index"
    hisat2_index_command = (
        f"~/TOOLS/hisat2-2.2.1/hisat2-build -p 32 {genome} hisat2.index"
    )
    check = cmd_run_multiple([hisat2_index_command], ou_path, check="hisat2_index")

    # run hisat2 alignment
    if check.exists():
        rna_lst = list(Path(rna_path).glob("*_1*fq.gz"))
        align_command_lst = []
        for PE1 in rna_lst:
            PE2 = Path(str(PE1).replace("_1.", "_2."))
            sam_out = ou_path / f"{str(PE1.name.split('_1.')[0])}.sam"
            hisat2_align = f"~/TOOLS/hisat2-2.2.1/hisat2 -x {hisat2_index} -1 {PE1} -2 {PE2} -p 32 --n-ceil L,0,0.05 --no-discordant --no-mixed --dta -S {sam_out}\n"
            align_command_lst.append(hisat2_align)

        # mutipleprocesses to run hisat2 alignment
        check = cmd_run_multiple(
            align_command_lst, ou_path, check="hisat2_align", num=8
        )
    else:
        print(f"CAN NOT FOUND CHECK `{check}` OF PREVIOUS STEP, SOMETHING HAS WRONG")
        exit(1)

    if check.exists():
        check = sam2bam(ou_path)
    else:
        print(f"CAN NOT FOUND CHECK `{check}` OF PREVIOUS STEP, SOMETHING HAS WRONG")
        exit(1)

def RNA_stringtie(ou_path, genome, verbose):
    root_path = Path(ou_path)
    # merging the bam of each tissues firstly
    ou_path = Path(root_path) / "RNA" / "hisat"
    check, merge_path = bam_merge(ou_path)

    if check.exists():
        # performing stringtie prediction
        ou_path = Path(root_path) / "RNA" / "stringtie"
        ou_path.mkdir(exist_ok=True, parents=True)
        stringtie_gtf = ou_path / "stringtie_prediction.gtf"
        stringtie_predict = f"~/TOOLS/stringtie-2.2.1.Linux_x86_64/stringtie -p 128 -o {stringtie_gtf} {merge_path}"
        check = cmd_run_multiple([stringtie_predict], ou_path, "stingtie_prediction")
    else:
        print(f"CAN NOT FOUND CHECK `{check}` OF PREVIOUS STEP, SOMETHING HAS WRONG")
        exit(1)

    # Finding coding regions using transdecoder
    ## The prepare work of running transdecoder
    if check.exists():
        ou_path = Path(root_path) / "RNA" / "TransDecoder"
        ou_path.mkdir(exist_ok=True, parents=True)
        cdna_path = ou_path / "transcripts.fasta"
        gff3_path = ou_path / "transcripts.gff3"

        cdna_command = f"/home/panda2bat/TOOLS/miniconda3/envs/transdecoder/opt/transdecoder/util/gtf_genome_to_cdna_fasta.pl {stringtie_gtf} {genome} > {cdna_path}"
        gff3_command = f"/home/panda2bat/TOOLS/miniconda3/envs/transdecoder/opt/transdecoder/util/gtf_to_alignment_gff3.pl {stringtie_gtf} > {gff3_path}"

        check = cmd_run_multiple(
            [cdna_command, gff3_command], ou_path, check="transdecoder_prepar"
        )
    else:
        print(f"CAN NOT FOUND CHECK `{check}` OF PREVIOUS STEP, SOMETHING HAS WRONG")
        exit(1)
    ## The Step 1: extract the long open reading frames
    if check.exists():
        step1_command = f"/home/panda2bat/TOOLS/miniconda3/envs/transdecoder/opt/transdecoder/TransDecoder.LongOrfs -t {cdna_path}"
        check = cmd_run_multiple(step1_command, ou_path, check="transdecoder_step1")
    else:
        print(f"CAN NOT FOUND CHECK `{check}` OF PREVIOUS STEP, SOMETHING HAS WRONG")
        exit(1)

    ## The Step 2: Including homology searches as ORF retention criteria
    if check.exists():
        uniref90 = Path("/home/panda2bat/DATABASE/UniRef90/uniref90.fasta")
        # makeblast_path = ou_path / "blastdb" / 'uniref90.db'
        # makeblast_path.mkdir(exist_ok=True, parents=True)
        # makeblastdb = f"~/TOOLS/ncbi-blast-2.11.0+/bin/makeblastdb -dbtype prot -in {uniref90} -out {makeblast_path}"

        # Instead blastp to diamond to speed up
        diamond_db_path = ou_path / "diamond_db" / "uniref90.db"
        diamond_db_path.mkdir(exist_ok=True, parents=True)
        diamond_db = f"~/TOOLS/diamond-2.0.14/bin/diamond makedb --in {uniref90} --db {diamond_db_path}"

        # check1 = cmd_run_multiple(makeblastdb, ou_path, check="makeblastdb", verbose=True)
        check1 = cmd_run_multiple(
            diamond_db, ou_path, check="diamond_makedb", verbose=verbose
        )

        if check1.exists():
            pfam_data = Path(
                "/home/panda2bat/Avivorous_bat/input/Download/Pfam/Pfam-A.hmm"
            )
            pfam_out = ou_path / "pfam.domtblout"
            pfam_search = f"/home/panda2bat/TOOLS/hmmer-3.3.2/bin/hmmsearch --cpu 128 -E 1e-10 --domtblout {pfam_out} {pfam_data} {ou_path}/longest_orfs.pep"

            # blastp = f"~/TOOLS/ncbi-blast-2.11.0+/bin/blastp -query {ou_path}/longest_orfs.pep -db {makeblast_path} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 72 > {ou_path}/blastp.outfmt6"
            # check = cmd_run_multiple(
            #    [pfam_search, blastp], ou_path, check="transdecoder_step2"
            # )
            # check = cmd_run_multiple(blastp, ou_path, check = "transdecoder_step2")

            diamand = f"~/TOOLS/diamond-2.0.14/bin/diamond blastp --db {diamond_db_path} --query {ou_path}/longest_orfs.pep --max-target-seqs 1 --outfmt 6 --threads 72 --sensitive --out {ou_path}/blastp.outfmt6 --evalue 1e-5"
            check = cmd_run_multiple(
                diamand, ou_path, check="transdecoder_step2", verbose=verbose
            )
        else:
            print("PLEASE perform `makeblastdb` firstly")
            exit(1)
    else:
        print(f"CAN NOT FOUND CHECK `{check}` OF PREVIOUS STEP, SOMETHING HAS WRONG")
        exit(1)

    ## The Step 3: predict the likely coding regions
    if check.exists():
        transdecoder_step3 = f"/home/panda2bat/TOOLS/miniconda3/envs/transdecoder/bin/TransDecoder.Predict -t {cdna_path} --retain_pfam_hits {pfam_out} --retain_blastp_hits {ou_path}/blastp.outfmt6"
        check = cmd_run_multiple(
            transdecoder_step3, ou_path, check="transdecoder_step3", verbose=verbose
        )
    else:
        print(f"CAN NOT FOUND CHECK `{check}` OF PREVIOUS STEP, SOMETHING HAS WRONG")
        exit(1)

    if check.exists():
        transdecoder_step4 = f"/home/panda2bat/TOOLS/miniconda3/envs/transdecoder/opt/transdecoder/util/cdna_alignment_orf_to_genome_orf.pl {ou_path}/transcripts.fasta.transdecoder.gff3 {ou_path}/transcripts.gff3 {ou_path}/transcripts.fasta > {ou_path}/transcripts.fasta.transdecoder.genome.gff3"
        cmd_run_multiple(
            transdecoder_step4, ou_path, check="transdecoder_step4", verbose=verbose
        )
    else:
        print(f"CAN NOT FOUND CHECK `{check}` OF PREVIOUS STEP, SOMETHING HAS WRONG")
        exit(1)

def trinity_pasa(ou_path, ss, verbose, genome, db_name):
    root_path = Path(ou_path)
    # Trinity genome guided assembly
    ou_path = root_path / "RNA" / "trinity_out"
    # print(ou_path)
    ou_path.mkdir(exist_ok=True, parents=True)
    BAM_path = root_path / "RNA" / "hisat" / "merge.bam"

    trinity_run = f"/home/panda2bat/TOOLS/miniconda3/envs/trinity/bin/Trinity --max_memory 200G --genome_guided_bam {BAM_path} --genome_guided_max_intron 10000 --CPU 150 --output {ou_path} --inchworm_cpu 150 --bflyHeapSpaceMax 10G --bflyCPU 150"
    if ss != "false":
        trinity_run = trinity_run + f" --SS_lib_type {ss}"
    check = cmd_run_multiple(trinity_run, ou_path, check="trintiy_run", verbose=verbose)

    # The PASA annotaion line
    if check.exists():
        trinity_gg_path = ou_path / "Trinity-GG.fasta"
        ou_path = root_path / "RNA" / "pasa"
        ou_path.mkdir(exist_ok=True, parents=True)

        # replace the default config file of pasa
        ## generating a randam name of SQLite path to avoid replicate
        
        if not db_name:
            db_name = "".join(random.sample(string.ascii_letters, 5))
        
        if list(ou_path.glob('pasa_align_*')):
            pasa_db_path = list(ou_path.glob('pasa_align_*'))[0]
        else:
            pasa_db_path = Path(ou_path / f'pasa_align_{db_name}')
        
        ## save the conf in output directory
        pasa_template = open(
            Path(__file__).absolute().parent / ".pasa.alignAssembly.Template.txt", "r"
        ).read()
        pasa_template = pasa_template.replace("$$SQLPATH$$", str(pasa_db_path))
        pasa_assembly_template = open(ou_path / "pasa.alignAssembly.conf", "w")
        pasa_assembly_template.writelines(pasa_template)
        pasa_assembly_template.close()

        pasa_align = f"singularity exec ~/TOOLS/pasapipeline.v2.5.2-devb.simg /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c {ou_path}/pasa.alignAssembly.conf -C -R -g {genome} --ALIGNERS blat,minimap2 -t {trinity_gg_path} --CPU 120 ##"

        if ss != "false":
            pasa_align = pasa_align.replace("##", "--transcribed_is_aligned_orient")
        else:
            pasa_align = pasa_align.replace("##", "")

        #print(pasa_align)
        check = cmd_run_multiple(pasa_align, ou_path, check="pasa_align", verbose=verbose)


def braker(genome, protein, ou_path, verbose):
    # Using Braker to perform ab initio annotation
    root_path = Path(ou_path)
    ou_path = root_path / "ab_initio" / "braker3"
    ou_path.mkdir(exist_ok=True, parents=True)

    bam_path = root_path / "RNA" / "hisat"
    bam_lst = [str(i) for i in list(bam_path.glob("*.sort.bam"))]

    # the soft masked genome
    braker_command = f"singularity exec /home/panda2bat/TOOLS/braker3.sif braker.pl --genome={genome} --prot_seq={protein} --bam={','.join(bam_lst)} --workingdir={ou_path} --gff3 --threads=96"
    cmd_run_multiple(braker_command, ou_path, check="branker_abintio", verbose=verbose)

def galba(genome, protein, ou_path, verbose):
    # User Galba to train augustus
    root_path = Path(ou_path)
    ou_path = root_path / "ab_initio" / "galba"
    out_path.mkdir(exist_ok=True, parents=True)

    galba_command = f"singularity exec ~/TOOLS/galba.sif galba.pl --genome={genome} --prot_seq={protein} --threads=96"
    cmd_run_multiple(galba_command, ou_path, check="galba_abintio", verbose=verbose)


def modified_miniprot(miniprot_anno, modified_miniprot):
    modified = open(modified_miniprot, 'w')
    with open (miniprot_anno, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('#'):
                if 'mRNA' in line:
                    n=1
                    modified.write(line + '\n')
                if 'CDS' in line:
                    tmp=line.split('\t')
                    id=re.findall(r'Parent=(.*?);', tmp[-1])[0]+'.'+str(n)
                    n+=1
                    tmp[-1]='ID='+id+';'+tmp[-1]
                    modified.write('\t'.join(tmp) + '\n')
    modified.close()

#def EVM(ou_path, genome, verbose, threads):
#    root_path = Path(ou_path)
#    ou_path = Path(ou_path) / "EVM"
#    ou_path.mkdir(exist_ok=True, parents=True)
#    evm_input = ou_path / "evm_input"
#    evm_input.mkdir(exist_ok=True, parents=True)
#
#    weight_file = Path(__file__).parent / ".EVM.weight.txt"
#
#    # generating the input files for EVM
#    braker_anno = root_path / "ab_initio" / "braker3" / "braker.gff3"
#    pasa_anno = root_path / "RNA" / "pasa" / "pasa_align.pasa_assemblies.gff3"
#    miniprot_anno = list((root_path / "miniprot").glob("*.gff"))[0]
#    transdecoder_annp = (
#        root_path
#        / "RNA"
#        / "TransDecoder"
#        / "transcripts.fasta.transdecoder.genome.gff3"
#    )
#
#    ### Before running EVM, there are some preparation
#    # merge gff3 of braker3 and stringtie+trandsdecoder
#    run = subprocess.Popen(
#        f"cat {braker_anno} {transdecoder_annp} > {evm_input}/gene_predictions.gff3",
#        shell=True,
#        stdout=subprocess.PIPE,
#        stderr=subprocess.STDOUT,
#        text=True,
#    )
#    run.wait()
#    shutil.copy(weight_file, (evm_input / '.EVM.weight.txt'))
#    shutil.copy(pasa_anno, (evm_input / 'transcript_alignments.gff3'))
#    # Modified the gff3 of miniprot to used in EVM
#    modified_miniprot(miniprot_anno, (evm_input / 'protein_alignments.gff3'))
#    
#    
#    # partition the EVM inputs
#    evm_partition_command = f"singularity exec ~/TOOLS/EVidenceModeler.v2.1.0.simg EVidenceModeler --genome {genome} --weights {weight_file} --gene_predictions {evm_input}/gene_predictions.gff3 --protein_alignments {evm_input}/protein_alignments.gff3 --transcript_alignments {evm_input}/transcript_alignments.gff3 --segmentSize 100000 --overlapSize 10000 --CPU {threads} --sample_id evm_result"
#    check = cmd_run_multiple(
#        evm_partition_command, ou_path, check="EVidenceModeler_run", verbose=verbose
#    )
#
#    if check.exists():
#        print("╰(○'◡'○)╮: ALL EVidenceModeler annotation COMMANDS HAVE COMPLETED")
def EVM(ou_path, genome, verbose, threads):
    # Do not use the results of pasa
    root_path = Path(ou_path)
    ou_path = Path(ou_path) / "EVM"
    ou_path.mkdir(exist_ok=True, parents=True)
    evm_input = ou_path / "evm_input"
    evm_input.mkdir(exist_ok=True, parents=True)

    # EVM weight 文件模板, 参考官网模板
    # The EVM.weight2.txt 移除了转录组证据
    weight_file = Path(__file__).absolute().parent / ".EVM.weight2.txt"

    # generating the input files for EVM
    braker_anno = root_path / "ab_initio" / "braker3" / "braker.gff3"
    pasa_anno = root_path / "RNA" / "pasa" / "pasa_align.pasa_assemblies.gff3"
    miniprot_anno = list((root_path / "miniprot").glob("*.gff"))[0]
    transdecoder_anno = (
        root_path
        / "RNA"
        / "TransDecoder"
        / "transcripts.fasta.transdecoder.genome.gff3"
    )

    ### Before running EVM, there are some preparation
    # merge gff3 of braker3 and stringtie+trandsdecoder
    run = subprocess.Popen(
        f"cat {braker_anno} {transdecoder_anno} > {evm_input}/gene_predictions.gff3",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    run.wait()
    shutil.copy(weight_file, (evm_input / '.EVM.weight2.txt'))
    # Modified the gff3 of miniprot to used in EVM
    modified_miniprot(miniprot_anno, (evm_input / 'protein_alignments.gff3'))
    
    
    # partition the EVM inputs
    evm_partition_command = f"singularity exec ~/TOOLS/EVidenceModeler.v2.1.0.simg EVidenceModeler --genome {genome} --weights {weight_file} --gene_predictions {evm_input}/gene_predictions.gff3 --protein_alignments {evm_input}/protein_alignments.gff3 --segmentSize 500000 --overlapSize 10000 --CPU {threads} --sample_id evm_result"
    check = cmd_run_multiple(
        evm_partition_command, ou_path, check="EVidenceModeler_run", verbose=verbose
    )

    if check.exists():
        print("╰(○'◡'○)╮: ALL EVidenceModeler annotation COMMANDS HAVE COMPLETED")

def pasa_update(ou_path, genome, verbose, threads):
    root_path = Path(ou_path)
    ou_path = root_path / 'PASA_update'
    ou_path.mkdir(exist_ok=True, parents=True)
    
    # Load the EVM annotation
    alignAssembly_path = root_path / 'RNA' / 'pasa' / 'pasa.alignAssembly.conf'
    evm_gff = root_path / 'EVM' / 'evm_result.EVM.gff3'
    Load_command = f"singularity exec ~/TOOLS/pasapipeline.v2.5.2-devb.simg /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi -c {alignAssembly_path} -g {genome} -P {evm_gff}"
    check = cmd_run_multiple(Load_command, ou_path, check = 'Load_Current_Gene_Annotations', verbose = verbose)
    
    # perfoming annotCompare update
    if check.exists():
        database_parameter = [i.strip() for i in open(alignAssembly_path).readlines() if i.startswith('DATABASE')][0]
        annotCompare_config = open(ou_path / 'pasa.annotCompare.conf', 'w')
        annotCompare_config.write(
            open(Path(__file__).absolute().parent / '.pasa.annotCompare.Template.txt', 'r').read().replace('DATABASE=$$', database_parameter)
        )
        annotCompare_config.close()
        transcript = root_path / 'RNA' / 'trinity_out' / 'Trinity-GG.fasta'
        
        annotCompare_command = f"singularity exec ~/TOOLS/pasapipeline.v2.5.2-devb.simg /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c {ou_path}/pasa.annotCompare.conf -A -g {genome} -t {transcript} --CPU {threads}"
        check = cmd_run_multiple(annotCompare_command, ou_path, check = 'annotCompare', verbose = verbose)
        
    if check.exists():
        print("╰(○'◡'○)╮: ALL PASA upsate COMMANDS HAVE COMPLETED") 

def longest_pep(gff):
    gene_data = {'parent':[]}
    mrna_data = {'parent':[], 'mrna':[]}
    prot_data = {'parent':[], 'mrna':[], 'seq':[], 'length':[]}
    longest_lst = []
    with open (gff, 'r') as file:
        for line in file:
            if not (line.startswith('#') or line == '\n'):
                line = line.strip()
                array = line.split('\t')
                if array[2] == 'gene':
                    gene_data['parent'].append(array[8].split(';')[0].replace('ID=',''))
                elif array[2] == 'mRNA':
                    mrna_data["parent"].append(array[8].split(';')[1].replace('Parent=',''))
                    mrna_data["mrna"].append(array[8].split(';')[0].replace('ID=', ''))
            elif line.startswith('#PROT'):
                line = line.strip()
                array = line.split('\t')
                prot_data['seq'].append(array[1])
                prot_data['length'].append(len(array[1]))
                array2 = array[0].split(' ')
                prot_data['parent'].append(array2[-1])
                prot_data['mrna'].append(array2[1])
    gene_df = pd.DataFrame(gene_data)
    mrna_df = pd.DataFrame(mrna_data)
    prot_df = pd.DataFrame(prot_data)
    
    info_df = gene_df.merge(mrna_df, on = 'parent', how = 'inner').merge(prot_df, on = ['parent', 'mrna'])
    info_df_max = info_df.groupby(by = 'parent').max()
    info_df_max = info_df_max.query('length >= 50')
    for index, data in info_df_max.iterrows():
        longest_lst.append(f'>{index}\n{data.seq.replace("*","")}\n')
    return(longest_lst)
  
def function_annotaion(ou_path, gff, verbose, threads):
    root_path = Path(ou_path)
    ou_path = root_path / 'Function'
    interpro_path = ou_path / 'InterPro'
    interpro_path.mkdir(exist_ok=True, parents=True)
    nr_path = ou_path / 'NR'
    nr_path.mkdir(exist_ok=True, parents=True)
    uniprot_path = ou_path / 'UniProt'
    uniprot_path.mkdir(exist_ok=True, parents=True)
    kofam_path = ou_path / 'KoFam'
    kofam_path.mkdir(exist_ok=True, parents=True)
    
    
    #obtain the longest sequence of gene
    if not (ou_path / 'longest_protein.fa').exists():
        longest_seq = longest_pep(gff)
        longest_seq_file = open(ou_path / 'longest_protein.fa', 'w')
        longest_seq_file.writelines(longest_seq)
        longest_seq_file.close()
    
    #perfome blast with NR and SwissProt database
    blast_NR_command = f"~/TOOLS/ncbi-blast-2.11.0+/bin/blastp -db /mnt/Extra_storage/NR/nr -query longest_protein.fa -evalue 1e-5 -outfmt 6 -num_alignments 1 -max_hsps 1 -num_threads 96 -out {nr_path}/blastp.nr.out"
    blast_UniProt_command = f"~/TOOLS/ncbi-blast-2.11.0+/bin/blastp -db /mnt/Extra_storage/UniProt/UniProt -query longest_protein.fa -evalue 1e-5 -outfmt 6 -num_alignments 1 -max_hsps 1 -num_threads {threads} -out {uniprot_path}/blastp.uniprot.out"
    #check = cmd_run_multiple([blast_UniProt_command], ou_path, check = 'blast_UniProt', verbose = verbose)
    print(blast_UniProt_command)
    
    if check.exists():
        #perfome InterProScan search
        interProScan_command = f"/mnt/Extra_storage/InterPro/interproscan-5.62-94.0/interproscan.sh -i longest_protein.fa -goterms --cpu {threads} -f tsv -o {interpro_path}/interProScan.tsv"
        check = cmd_run_multiple(interProScan_command, ou_path, check = 'InterProScan', verbose = verbose)
    
    if check.exists():
        #perfome Kofam KEGG pathways search
        kofam_command = f"/home/panda2bat/TOOLS/kofam_scan-1.3.0/bin/exec_annotation -c /home/panda2bat/TOOLS/kofam_scan-1.3.0/bin/config.yml -E 0.00001 --tmp {kofam_path}/tmp -f mapper -o {kofam_path}/kofam.results longest_protein.fa"
        check = cmd_run_multiple(kofam_command, ou_path, check = 'KoFam', verbose = verbose)
    
def main():
    args = make_parse()
    step = args.step
    genome = args.genome
    ou_path = args.ou_path
    verbose = args.verbose

    if step == 1:
        prefix = args.prefix
        protein = args.protein
        miniprot_run(genome, ou_path, protein, prefix, verbose)
    elif step == 2:
        rna_path = args.rna
        RNA_anno(rna_path, ou_path, genome, verbose)
    elif step == 3:
        RNA_stringtie(ou_path, genome, verbose)
    elif step == 4:
        ss = args.strand
        dbname = args.dbname
        trinity_pasa(ou_path, ss, verbose, genome, dbname)
    elif step == 5:
        protein = args.protein
        braker(genome, protein, ou_path, verbose)
        #galba(genome, protein, ou_path, verbose)
    elif step == 6:
        threads = args.threads
        EVM(ou_path, genome, verbose, threads)
    elif step == 7:
        threads = args.threads
        pasa_update(ou_path, genome, verbose, threads)
    elif step == 8:
        threads = args.threads
        gff = args.gff
        function_annotaion(ou_path, gff, verbose, threads)
    else:
        print("Please the correct step code")
        exit(1)


if __name__ == "__main__":
    main()
    
