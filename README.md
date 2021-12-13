# RMetaflow (RMF) - An R-Based Pipeline for Metagenomic Processing 

## Introduction 
RMetaflow uses R scripts to guide sequence processing and analysis coordinated with the use of singularity containers.

*Why another metagenomics pipeline?* Other pipelines exist but lack some control and flexibility that I was seeking for my own analysis. Other pipelines have great merits and may be better for you.

*Why R?* R efficiently oversees the pipeline, is easy to code, is readable, and has many tools for downstream data analysis, alleviating some need to jump back and forth between programming languages.  Since the heavy computational steps are run in singularity containers, R is not itself using major resources. 

The major steps for RMetaflow include:
1. QC: sequence filtering, trimming, decontamination, and deduplication. 
   
   PROGRAMS: 
   1. fastp
   2. CONCOCT cut_up_fasta.py 
   3. minimap2 
   4. samtools 
   5. bbtools 
   6. pigz
   7. idseq-dedup 
   
2. Classification of paired reads. 
   
   PROGRAMS: 
   1. kaiju
   2. humann (if you like)
   
3. De novo assembly of the reads, classification, and functional analysis.
   
    PROGRAMS: 
    1. megahit  
    2. prinseq-lite 
    3. metaerg 
    4. CoverM 
    5. minpath 
   
## Software Requirements
1. A singularity build environment (although you may use faker root if necessary)
2. R. Expected to work on version >= 3.6  but routinely run using version 4.1
3. Text editing software (I do this right in RStudio for convenience but not necessary)
4. Build singularity container qc.sif from script qc.def. Place the sif file in the containers folder.
5. Build singularity container classify.sif from script classify.def. Place the sif file in the containers folder.
6. Build singularity container assembly-by-refseq.sif from script assemby-by-refseq.def. Place the sif file in the containers folder.
7. signalp-5.0b.Linux.tar.gz (in the build directory) from: https://services.healthtech.dtu.dk/cgi-bin/sw_request
8. FragGeneScan1.31.tar.gz (in the build directory) from: https://sourceforge.net/projects/fraggenescan/files/
9. MaxBin-2.2.7.tar.gz (in the build directory) from: https://sourceforge.net/projects/maxbin/files/

## Database Requirements
1. For Kaiju: 
https://kaiju.binf.ku.dk/database/kaiju_db_progenomes_2020-05-25.tgz
* Place it in databases/kaiju  
* other databases are also available including all of nr. Your run time for classification will increase with larger databases 
2. For decontamination:
* make directories and download references
  * databases/contaminants/human 
    * wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
  * databases/contaminants/mouse  
    * wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
  * databases/contaminants/phiX 
    * wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz
  * gunzip -d each of the files in the respective folders
3. For Humann use the **classify.sif** container to run the following:
* humann_databases --download chocophlan full databases/humann --update-config yes
* humann_databases --download uniref uniref90_diamond /databases/humann --update-config yes
* humann_databases --download utility_mapping full /databases/humann --update-config yes
4. For Metaerg:
* databases/metaerg 
  * wget -c http://ebg.ucalgary.ca/metaerg/db.tar.gz
  * tar xvfz db.tar.gz 

## The Yaml File 
1. You will need a folder "yaml" containing the file "metaflow-steps.yml" that will guide each of the processing steps.
2. Each of the steps has some basic components 
   1. The step name which is a combination of the module and step (module.step). Note the position of the period and end colon 
   2. from: this is the data source for that step
   3. in_ext: the file extension of the input files 
   4. out_ext: the expected output file extension 
   5. pass: usually empty but there to pass additional variables
   6. source: the script with the step commands. Note that the source has 2 parts.  First is either "source" or "slurm" (if running with a slurm server). These are separated by 3 x ":" like this source:::R/script.R 
   7. completed: this is 'yes' or 'no' to indicate if the step was completed. You can manually change this is in the yaml file if you want to skip steps or rerun them.
3. The first step must be "startup:" This step has the unique field of "split_ext" to tell the program where to split the files names and obtain the leftmost part of the split to designate as sample names.  There can be multiple splits using "|" as a divider indicating 'or'.  Like this "_L1|_L2" will split the file at an "_L1" or "_L2" position. If the file is sample01_L1_R1.fastq.gz, the sample name would be sample01. 

## Run conditions
1. I routinely run the pipeline and containers on a Slurm-based management system where R and singularity are loaded as modules. The use of slurm can be controlled. 
2. I routinely will run with arrays of 100 nodes and up to 52 threads per node. Ideally you have at least 64 Gb per node available but 128 Gb per node is better. The number of nodes and threads can be modified.

## Citations
* Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018 Sep 1;34(17):i884–i890. PMCID: PMC6129281
* Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018 Sep 15;34(18):3094–3100. PMCID: PMC6137996
* idseq-dedup: https://github.com/chanzuckerberg/idseq-dedup
* Alneberg J, Bjarnason BS, de Bruijn I, Schirmer M, Quick J, Ijaz UZ, Lahti L, Loman NJ, Andersson AF, Quince C. Binning metagenomic contigs by coverage and composition. Nat Methods. 2014 Nov;11(11):1144–1146. PMID: 25218180
* Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078–2079. PMID: 19505943
* BBtools: BBMap – Bushnell B. – sourceforge.net/projects/bbmap/
* pigz:  https://zlib.net/pigz/
* Menzel P, Ng KL, Krogh A. Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nat Commun. 2016 Apr 13;7:11257. PMCID: PMC4833860
* Beghini F, McIver LJ, Blanco-Míguez A, Dubois L, Asnicar F, Maharjan S, Mailyan A, Manghi P, Scholz M, Thomas AM, Valles-Colomer M, Weingart G, Zhang Y, Zolfo M, Huttenhower C, Franzosa EA, Segata N. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. Elife [Internet]. 2021 May 4;10. Available from: http://dx.doi.org/10.7554/eLife.65088 PMCID: PMC8096432
* Li D, Liu C-M, Luo R, Sadakane K, Lam T-W. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics. 2015 May 15;31(10):1674–1676. PMID: 25609793
* Schmieder R, Edwards R. Quality control and preprocessing of metagenomic datasets. Bioinformatics. 2011 Mar 15;27(6):863–864. PMCID: PMC3051327
* Dong X, Strous M. An Integrated Pipeline for Annotation and Visualization of Metagenomic Contigs. Front Genet. 2019 Oct 15;10:999. PMCID: PMC6803454
* Woodcroft BJ. CoverM: Read coverage calculator for metagenomics [Internet]. Github; [cited 2021 Dec 13]. Available from: https://github.com/wwood/CoverM
* Wu Y-W, Simmons BA, Singer SW. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics. 2016 Feb 15;32(4):605–607. PMID: 26515820
