#' RMetagenomics
#' Initiated March 2021
#' Patrick C. Seed
#' Stanley Manne Children's Research Institute
#' Ann and Robert H. Lurie Children's Hospital
#' Northwestern University
#' MIT license
#' 
#' Define variables

library(tidyverse)

# base machine information ------------------------------------------------

cores = threads = benchmarkme::get_cpu()
mem = benchmarkme::get_ram()
memory = round(mem * 0.8 /10E8, 0)

printGlue("#----------------------------------------------------------------")
printGlue("#----Running RMetaflow with ")
printGlue("#----{cores} cores available ")
printGlue("#----{memory}GB = 80% of memory available ")
printGlue("#----------------------------------------------------------------")

# sample and input information --------------------------------------------

projectdir = getProjDir()
basedir = "/mnt" # basedir = dirname(projectdir)
projectname = basename(projectdir)
  
current.step = commands[[1]]
in_ext <<- current.step$in_ext
split_ext <<- current.step$split_ext

if (current.step$pair_ext == "_1"){
  pair_ext <<- c("_1", "_2")
}
if (current.step$pair_ext == ".R"){
  pair_ext <<- c(".R1", ".R2")
}else{
  pair_ext <<- c("_R1", "_R2")
}

samples <<- sapply(list.files(dirNormalize(current.step$from), pattern = current.step$in_ext), function(x) str_split(x, split_ext)[[1]][1]) %>% as.vector() %>% unique()

# SLURM variables ---------------------------------------------------------

partition = "genomics"; account = "b1042"
#partition = "buyin"; account = "b1051"
mem = "96G"
time = "08:00:00" #this may change based on steps
N = 1
ntasks_per_node = threads = 28 #only allow 28 for b1051/buyin; can have up to 52 for b1042/genomics

# containers  --------------------------------------------------------------

singularity = stringr::str_glue("singularity run -B {basedir} ")
qc = stringr::str_glue("{singularity} containers/qc.sif ")
classify = stringr::str_glue("{singularity} containers/classify.sif ")
assembly.by.refseq = stringr::str_glue("{singularity} containers/assembly-by-refseq.sif ")
metaerg = stringr::str_glue("{singularity}containers/assembly-by-refseq.sif ")
coverm = stringr::str_glue("{singularity} containers/assembly-by-refseq.sif ")
semibin = stringr::str_glue("{singularity} run -B {basedir} containers/assembly-by-refseq.sif ")

fastp = stringr::str_glue("{singularity} containers/fastp.sif ")
bowtie = stringr::str_glue("{singularity} containers/bowtie2.sif ")
samtools = stringr::str_glue("{singularity} containers/samtools.sif ")
kaiju = stringr::str_glue("{singularity} containers/kaiju.sif ")
skesa = stringr::str_glue("{singularity} containers/skesa.sif ")
metaprokka = stringr::str_glue("{singularity} containers/metaprokka.sif ")
abricate = stringr::str_glue("{singularity} containers/abricate.sif ")
seqtk = stringr::str_glue("{singularity} containers/seqtk.sif ")
metahipmer = stringr::str_glue("{singularity} containers/metahipmer.sif ")
bbtools = stringr::str_glue("{singularity} containers/bbtools.sif ")
eggnogmapper = stringr::str_glue("{singularity} containers/eggnog-mapper.sif ")
squeezemeta = stringr::str_glue("{singularity} containers/squeezemeta.sif ")
pigz = stringr::str_glue("{singularity} containers/pigz.sif ")
idseq.dedup = stringr::str_glue("{singularity} containers/idseq-dedup.sif ")
decontaminate = stringr::str_glue("{singularity} containers/decontaminate.sif ")
minimap2 = stringr::str_glue("{singularity} containers/minimap2_v2.15dfsg-1-deb_cv1.sif")
seqstats = stringr::str_glue("{singularity} containers/seqstats.sif ")
flash = stringr::str_glue("{singularity} containers/pangenome.sif ")
mash = stringr::str_glue("{singularity} containers/pangenome.sif ")
sambamba = stringr::str_glue("{singularity} containers/pangenome.sif ")
cdhit = stringr::str_glue("{singularity} containers/squeezemeta.sif ")
megahit = stringr::str_glue("{singularity} containers/pangenome.sif ")
orna = stringr::str_glue("{singularity} containers/orna.sif ")
quickmerge = stringr::str_glue("{singularity} containers/pangenome.sif ")
prinseq = stringr::str_glue("{singularity} containers/squeezemeta.sif ")
humann = stringr::str_glue("{singularity} containers/humann3.sif ")
seqkit = stringr::str_glue("{singularity} containers/seqkit.sif")
MAC = stringr::str_glue("{singularity} containers/MAC.sif")
flye = stringr::str_glue("{singularity} containers/flye.sif ")
abricate = stringr::str_glue("{singularity} containers/abricate.sif ")

# MODULE-SPECIFIC VARIABLES -----------------------------------------------

variables.fastp = stringr::str_glue(" -p -w 16 ")
variables.semibin = stringr::str_glue(" -m 1000 --ml-threshold  --random-seed 1234 ")
variables.metaerg = stringr::str_glue(" --cpus {threads} --dbdir databases/metaerg/db  ")
variables.bowtie.decon = stringr::str_glue("-x databases/contaminants --threads {threads} --very-sensitive-local ")
variables.kaiju.map = stringr::str_glue(" -z {threads} -f databases/kaiju/progenomes/kaiju_db_progenomes.fmi -t databases/kaiju/progenomes/nodes.dmp -a mem")
variables.kaiju.tax = stringr::str_glue( " -n databases/kaiju/progenomes/names.dmp -t databases/kaiju/progenomes/nodes.dmp -p ")
variables.kaiju.table = stringr::str_glue(" -n databases/kaiju/progenomes/names.dmp -t databases/kaiju/progenomes/nodes.dmp -r species -l superkingdom,phylum,class,order,family,genus,species -m 0.001 ")
variables.kaiju.krona = stringr::str_glue(" -n databases/kaiju/progenomes/names.dmp -t databases/kaiju/progenomes/nodes.dmp ")
variables.skesa = stringr::str_glue(" --cores {cores} --memory {memory} ")
variables.abricate = stringr::str_glue(" ")
variables.idseq.dedup = stringr::str_glue("-l 70")
variables.metaprokka.merged = stringr::str_glue(" --cpus {threads} --mincontiglen 500 --force ")
variables.decontaminate = c("fna$", "databases/contaminants/human/", "databases/contaminants/mouse/",
                            "databases/contaminants/phiX/")
variables.squeezemeta.merged = stringr::str_glue(" -m merged -t {threads} --nopfam --extdb databases/squeezemeta ")
variables.squeezemeta.sequential = stringr::str_glue(" -m sequential -t {threads} --nopfam --extdb databases/squeezemeta ")
variables.pigz = stringr::str_glue(" -p {threads} ")
variables.megahit = stringr::str_glue(" -t {threads}  --memory 0.8 ")
variables.flash = stringr::str_glue(" -t 4 -M 150 ")
variables.orna = stringr::str_glue("-kmer 27 -base 1.7 -nb-cores {threads}")
variables.humann = stringr::str_glue("--threads {threads} --remove-temp-output --nucleotide-database databases/humann/chocophlan --protein-database databases/humann/uniref --metaphlan-options '--bowtie2db /opt/miniconda-4.6.14/envs/biobakery/lib/python3.7/site-packages/metaphlan/metaphlan_databases' ")
variables.flye = stringr::str_glue("--threads {threads} --iterations 0 --meta --subassemblies ")
variables.abricate = stringr::str_glue("--threads {threads} --db ncbi --nopath ")

methods = list(
  fastp =  "Sequences were quality controled and preprocessed using the software fastp with the default settigs {references$fastp}.",
  bbduk = "Human and mouse sequences were filtered out using the tool bbduk using the references GRCh38.p13 and GRCm39, respectively {references$bbduk}.",
  bowtie = "Human, mouse, and phiX174 sequences were filtered out using bowtie2 and the sequence references GRCh38.p13, GRCm39, and GCF_000819615.1, respectively {references$bowtie}.",
  flash = "Pre-processed paired end reads were merged by overlap using the tool FLASH {references$flash}.",
  mmseqs= "Taxonomy was assigned by the 2bLCA method using the tool mmseq2 and the reference database {references$mmseq}.",
  kaiju = "Taxonomy was assigned using the tool kaiju and the reference database {references$kaiju}.",
  megahit = "Pre-processed metagenomic sequences were de novo assembled using the tool MEGAHIT {references$megahit}.",
  kneaddata = "",
  skesa = " Pre-processed metagenomic sequences were de novo assembled using the tool SKESA {references$skesa}. ",
  abricate = "", 
  seqtk = "",
  metahipmer = "",
  bbtools = "",
  squeezemeta = "",
  pigz = ""
)

# DATABASE REFERENCES -----------------------------------------------------

database.references = list()

# WHERE TO PULL CONTAINERS ------------------------------------------------

containers.pull = list(
  samtools =list("biocontainers/samtools:v1.9-4-deb_cv1" , "samtools.sif"),
  fastp = list( "biocontainers/fastp:v0.20.1_cv1", "fastp.sif" ),
  bowtie2 =  list("biocontainers/bowtie2:v2.4.1_cv1", "bowtie2.sif" ),
  kaiju = list("nanozoo/kaiju:1.7.2--fa366a0", "kaiju.sif"),
  skesa = list( "staphb/skesa:2.4.0", "skesa.sif"),
  abricate = list("staphb/abricate:1.0.0", "abricate.sif"),
  seqk = list("biocontainers/seqtk:v1.3-1-deb_cv1", "seqtk.sif"),
  bbtools = list("bryce911/bbtools:38.90", "bbtools.sif"),
  pigz = list("rtibiocloud/pigz:v2.4_b243f9", "pigz.sif")
)

containers.build = list(
  metaprokka = list("definitions/metaprokka.def", "metaprokka.sif"),
  r_base = list("definitions/r-base.def", "r-base.sif"),
  eggnog_mapper = list("definitions/eggnog-mapper.def", "eggnog-mapper.sif"),
  squeezemeta = list("definitions/squeezemeta.v2.def", "squeezemeta.sif"),
  seqstats = list("definitions/seqstats.def", "seqstats.sif")
)

# decorations
indent1 = "#---- "
indent2 = "#-------- "
indent3 = "#------------ "
indent4 = "#---------------- "

