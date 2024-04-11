# A chromosome-level genome assembly of a avivorous bat (Nyctalus aviator

Currently, there are three carnivorous bat species, namely *Ia io*, *Nyctalus lasiopterus*, and *Nyctalus aviator*, are known to actively prey on seasonal migratory birds (hereinafter referred to as “avivorous bats”). However, the absence of reference genomes impedes a thorough comprehension of the molecular adaptations of avivorous bat species. Herein, we present the high-quality chromosome-scale reference genome of *N. aviator* based on PacBio subread, DNBSeq short-read and Hi-C sequenceing data. The genome assembly size of *N. aviator* is 1.77 Gb, with a scaffold N50 of 102 Mb, of which 99.8% assembly was  anchored into 20 pseudo-chromosomes. After masking 514.3 Mb repetitive sequences, a total of 19,703 protein-coding genes were identified, of which 96.2% were functionally annotated. The genome assembly and gene prediction reached 95.6% and 96.5% completeness of Benchmarking Universal Single-Copy Orthologs (BUSCO), respectively. This chromosome-level reference genome of *N. aviator* fills a gap in the existing information on the genomes of carnivorous bats, especially avivorous ones, and will be valuable for mechanism of adaptations to dietary niche expansion in bat species.

Detailed information about the parameters and custom scripts utilized in this artical can be obtained from this repository 

## Code/Software

`Python >= 3.8.10`

Dependence:
`pandas`

`Biopython`

`rpy2`

`tqdm`

`gffutils`

Software:
[`fastp`](https://github.com/OpenGene/fastp) 

[`jellyfish`](https://github.com/gmarcais/Jellyfish)

[`nextDenovo`](https://github.com/Nextomics/NextDenovo)

[`wtdbg2`](https://github.com/ruanjue/wtdbg2)

[`nextPolish`](https://github.com/Nextomics/NextPolish)

[`YaHS`](https://github.com/c-zhou/yahs)

[`busco`](https://busco.ezlab.org)

[`EDTA`](`07.genome-repeat.py`)

[`hisat2`](https://daehwankimlab.github.io/hisat2/)

[`stringtie`](https://ccb.jhu.edu/software/stringtie/)

[`TransDecoder`](https://github.com/TransDecoder/TransDecoder/wiki)

[`PASA`](https://github.com/PASApipeline/PASApipeline)

[`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

[`braker3`](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

[`miniprot`](https://github.com/lh3/miniprot)

[`EVM`](https://github.com/EVidenceModeler/EVidenceModeler/wiki#alternatives-to-evm)

[`LASTZ`](https://www.bx.psu.edu/~rsharris/lastz/)

[`CNEr`](https://bioconductor.org/packages/release/bioc/html/CNEr.html)

[`make_lastz_chains`](https://github.com/hillerlab/make_lastz_chains)

[`inparanoid`](https://inparanoidb.sbc.su.se)

[`iqtree2`](http://www.iqtree.org)

[`ASTRALIII`](https://github.com/smirarab/ASTRAL)

[`pal2nal`](https://www.bork.embl.de/pal2nal/)

[`prank`](https://github.com/ariloytynoja/prank-msa)

[`samtools`](https://www.htslib.org)

**Prior to usage, verify the proper installation of all dependencies and software. As the scripts underwent testing solely in a local environment, adjusting the software paths within the script may be required based on the user's environment.**

## Usage
### 1. Quality control
The script `01.quality.py` utilizes [`fastp`](https://github.com/OpenGene/fastp) to filter  sequences containing over 20% low-quality bases or more than 10% 'N' bases.  
```
usage: 01.quality.py [-h] --flist FLIST [--in INPUT] [--ou OUPUT] [--na NAME]

optional arguments:
  -h, --help     show this help message and exit
  --flist FLIST  The file contains names of all input files. The file contains two columns, the delimter is tab. The first column is the name of first file of PE paried, the second is the second of PE paired. The file should not contain the paths of any input files, the path of input file can be set uing parameter `--in`
  --in INPUT     The path of directory, which contains all input files
  --ou OUPUT     The ouput path, If the output path is not set, it defaults to the input path.
  --na NAME      The prefix of output files
```

### 2. Genome Survey Analysis
The script `02.genome-survey.py` estimates the genome size from Illumina short reads using [`jellyfish`](https://github.com/gmarcais/Jellyfish) with a progressive k-mer approach.
```
usage: 02.genome-survey.py [-h] [--flist FLIST] [--in INPUT] [--ou OUTPUT] [--na NAME] [--range KRANGE] [--jpath JPATH]

optional arguments:
  -h, --help      show this help message and exit
  --flist FLIST   The file contains all input, each row was divided two columns by 'tab', the first column was the first file of a paired-end fastq file, and the second row is the second
  --in INPUT      the path of the directory, which contains all fastq files
  --ou OUTPUT     the path of output directory of resutls
  --na NAME       the prefix of results files
  --range KRANGE  The program iterates through a range of k-mer sizes specified as start:end:step, such as 10:21:1 The default k-mer range spans from 10 to 21 with a step size of 1.
  --jpath JPATH   The absoult path of jellyfish programe
```

### 3. Assembly
The script `03.genome-assembly.py` assembles PacBio subreads into contig-level assemblies using [`nextDenovo`](https://github.com/Nextomics/NextDenovo) and [`wtdbg2`](https://github.com/ruanjue/wtdbg2). This script follows this process: 1) acquiring consensus reads with `nextDenovo`; 2) obtaining paired alignments with `kbm2` in `wtdbg2`; 3) executing `wtdbg2` with paired alignments using  progressive parameters; 4) selecting the output with the highest contigN50 and genome size values as the optimal result.
```
usage: 03.genome-assembly.py [-h] [--in INPUT] [--ou OUTPUT] [--gsize GSIZE] [--step STEP]

optional arguments:
  -h, --help     show this help message and exit
  --in INPUT     the path of input directory, which contains all input files
  --ou OUTPUT    the path of output directory
  --gsize GSIZE  the prediceted size of genome, a string ended with g, for example, 2g
  --step STEP    the different step in process, 1: nextDenovo; 2: kbm2 in wtdbg2; 3: wtdbg in wtdbg2; 4 produce the best assembly
```
Users can modify the parameters of the wtdbg2 assembly by editing the script file (line: 113-116).
```
110     '''                                   
111     Modified parameters used in wtdbg2 assembly process 
112     ''' 
113     p_nodelen = [1536, 2048, 2304, 2560, 1024] 
114     p_s = [0.5] 
115     p_alndovetail = [-1, 4608, 9216, 256]  
116     p_L = [5000, 0]  
```


The script `04.genome-polish.py` performs polishing on the contig-level assembly derived from Illumina short reads and PacBio subreads using [`nextPolish`](https://github.com/Nextomics/NextPolish).
```
usage: 04.genome-polish.py [-h] [--initial INITIAL] [--inl LGS_PATH] [--ins SGS_PATH] [--ou OU_PATH]

optional arguments:
  -h, --help         show this help message and exit
  --initial INITIAL  the initial assembly fasta from wtdbg2
  --inl LGS_PATH     the input path for directory, which contains PacBio subreads
  --ins SGS_PATH     the input path for directory, which contains Illumina short reads
  --ou OU_PATH       the output directory of nextPolish
```

The script `05.genome-scaffold.py` scaffolds the contig-level assembly into scaffold-level using [`YaHS`](https://github.com/c-zhou/yahs). ~~In this script, you can select which tools be used, [`3D-DNA`](https://github.com/aidenlab/3d-dna) or [`YaHS`](https://github.com/c-zhou/yahs)~~
```
usage: 05.genome-scaffold.py [-h] [--genome GENOME] [--hic_path HICPATH] [--ou_path OUPATH] [--step STEP] [--review REVIEW]

optional arguments:
  -h, --help          show this help message and exit
  --genome GENOME     The path of polished fasta file
  --hic_path HICPATH  The path of directory, which contains all hic fastq files
  --ou_path OUPATH    The path of output directory
  --step STEP         Each number corresponds to a specific step, 1:chromap index; 2:chromap alignment; 3:yash scaffold; 4:yash2juicer; 5:manual curation by juicer; 6: review; 7: generate the hic concat map file
  --review REVIEW     The review results of juicer, only be used in second round of 3D-DNA
```

The script `06.genome-completeness.py` evaluates genome completeness utilizing [`busco`](https://busco.ezlab.org/).
```
usage: 06.genome-completeness.py [-h] [--genome GENOME] [--ou OU_PATH] [--protein PROTEIN] [--step STEP]

optional arguments:
  -h, --help         show this help message and exit
  --genome GENOME    The finally assembly genome fasta file
  --ou OU_PATH       The path of output directory
  --protein PROTEIN  The protein fasta of annotation, only be used in annotaion completes analysis
  --step STEP        Each number corresponds to a specific step: 1:busco with metaeuk; 2: busco with augustus; 3: annotation completes using busco
```

### 4. Annotation

The script `07.genome-repeat.py` annotates repeat sequences by using [`EDTA`](`07.genome-repeat.py`).
```
usage: 07.genome-repeat.py [-h] [--genome GENOME] [-o OU_PATH] [--cds CDS] [--curatedlib CURATED] [--step STEP]

optional arguments:
  -h, --help            show this help message and exit
  --genome GENOME       The finally genome assembly
  -o OU_PATH            The path of output directory
  --cds CDS             The cds sequences of closely-related specie
  --curatedlib CURATED  The curated lib file
  --step STEP           Each number corresponds to a specific step, 1:EDTA; 2:softmask;
```

Various annotation tools were utilized in `08.genome-annotation.py` to annotate the genome structures of the assembly. The tools employed include [`hisat2`](https://daehwankimlab.github.io/hisat2/), [`stringtie`](https://ccb.jhu.edu/software/stringtie/), [`TransDecoder`](https://github.com/TransDecoder/TransDecoder/wiki), [`PASA`](https://github.com/PASApipeline/PASApipeline), [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki), [`braker3`](https://github.com/trinityrnaseq/trinityrnaseq/wiki), [`miniprot`](https://github.com/lh3/miniprot), and [`EVM`](https://github.com/EVidenceModeler/EVidenceModeler/wiki#alternatives-to-evm).
```
usage: 08.genome-annotation.py [-h] [-g GENOME] [-p PROTEIN] [-o OU_PATH] [-n PREFIX] [--rna RNA] [--SS {false,RF,FR}] [--dbname DBNAME] [--gff GFF] [-s STEP] [--verbose VERBOSE] [--threads THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        The path of genome input (softmasked)
  -p PROTEIN, --protein PROTEIN
                        The path of protein fasta file of closely-related species, this protein file will be used in miniprot and braker3
  -o OU_PATH, --ou_path OU_PATH
                        The path of output directory
  -n PREFIX, --name PREFIX
                        The prefix of miniprot output file
  --rna RNA             The path of directory, which contains RNA fastq files
  --SS {false,RF,FR}    The BOOL value of `SS_lib_type` parameter in Trintiy
  --dbname DBNAME       The SQL name used in pasa, if not gived, the pasa db will named with 5 random letters
  --gff GFF             The GFF annotation of PASA update, this option only used in step 8.
  -s STEP, --step STEP  The step in annotation process. 1:miniprot; 2:hisat2; 3:stringtie; 4:trinity+pasa; 5:braker; 6:EvidenceModeler; 7:PASA upate; 8:functional annotation
  --verbose VERBOSE     Print the verbpse information (True) or not (Flase), default is False
  --threads THREADS     The number of threads
```

### 5. Systeny
The script `09.genome-systeny.py` generates pairwise genome alignments with [`LASTZ`](https://www.bx.psu.edu/~rsharris/lastz/). To expedite the alignment, the `target` and `query` genome files were partitioned based on scaffolds using the R package [`CNEr`](https://bioconductor.org/packages/release/bioc/html/CNEr.html).
> It is recommended to use [`make_lastz_chains`](https://github.com/hillerlab/make_lastz_chains) to accelerate computation.

```
usage: 09.genome-systeny.py [-h] [--target TARGET] [--query QUERY] [--output OU_PATH] [--tprefix TPREFIX] [--qprefix QPREFIX] [--step STEP] [--verbose {True,False}] [--threads THREADS]

optional arguments:
  -h, --help            show this help message and exit
  --target TARGET, -t TARGET
                        The softmasked target genome used to perform lastz alignment. The genome files should be 2bit format
  --query QUERY, -q QUERY
                        The softmasked query genome, 2bit format
  --output OU_PATH, -o OU_PATH
                        The path of output directory
  --tprefix TPREFIX     The name of target
  --qprefix QPREFIX     The name of query
  --verbose {True,False}
                        Print all (True, default) information of programe or not (False)
  --threads THREADS     The nubmer of threads
```

### 6. Evolution
The script `10.evolution-obtain_the_longest_sequence_from_NCBI.py` is designed to retrieve the longest isoform for every coding gene. The script is capable of trimming the CDS sequences extracted from the gff3 file based on correspond protein sequence to ensure that all CDS sequence lengths are multiples of three. CDS sequences with a length shorter than three times the corresponding protein sequence will be removed.
```
usage: 10.evolution-obtain_the_longest_sequence_from_NCBI.py [-h] [-g GFF] [-f FNA] [--prefix PREFIX]

optional arguments:
  -h, --help         show this help message and exit
  -g GFF, --gff GFF  The NCBI gff3 file
  -f FNA, --fna FNA  The geneome file
  --prefix PREFIX    The prefix of the longest files, such as: prefix.gff3, prefix.cds, prefix.pep
```

The script `11.evolution-single_copy_gene_inparanoid.py` is designed to identified single-copy orthologous genes by using [`inparanoid`](https://inparanoidb.sbc.su.se)
```
usage: 11.evolution-single_copy_gene_inparanoid.py [-h] [-i INPUT] [-o OUTPUT] [-m MRNA] [--step STEP]

optional arguments:
  -h, --help   show this help message and exit
  -i INPUT     The input directory
  -o OUTPUT    The output directory
  -m MRNA      The corresponding mRNA directory
  --step STEP  Each number corresponds to a specific step: 1:inparanoid; 2: parse the results of inparanoid
```

The script `12.evolution-species_tree.py` is designed to reconstruct phylogenetic tree (species tree) based on the process: 1) MSA using [`prank`](https://github.com/ariloytynoja/prank-msa); 2) [`pal2nal`](https://www.bork.embl.de/pal2nal/); 3) [`iqtree2`](http://www.iqtree.org) on every orthologous gene; 4) [`ASTRALIII`](https://github.com/smirarab/ASTRAL) ; 5) estimated branch length with supergene matrix (concatenated orthologous genes) by using [`iqtree2`](http://www.iqtree.org) 
```
usage: 12.evolution-species_tree.py [-h] [-i INPUT] [-o OUTPUT] [--threads THREADS] [--step STEP] [-m MRNA] [-t TREE]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        The input directory
  -o OUTPUT, --output OUTPUT
                        The output directory
  --threads THREADS     The number of threads
  --step STEP           Each number corresponds to a specific step, 1) prank MSA; 2) pal2nal; 3) iqtree with every gene; 4) ASTRALIII; 5) estimated branch length with concatenated genes
  -m MRNA               The directory contains corresponding mRNA fasta files
  -t TREE               The labeled tree file used in mcmctree
```

