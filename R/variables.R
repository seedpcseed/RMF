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
mem = "96G"
time = "08:00:00" #this may change based on steps
N = 1
ntasks_per_node = threads = 28

# containers  --------------------------------------------------------------

singularity = stringr::str_glue("singularity run -B {basedir} ")
qc = stringr::str_glue("{singularity} containers/qc.sif ")
classify = stringr::str_glue("{singularity} containers/classify.sif ")
assembly.by.refseq = stringr::str_glue("{singularity} containers/assembly-by-refseq.sif ")
metaerg = stringr::str_glue("{singularity}containers/assembly-by-refseq.sif ")
coverm = stringr::str_glue("{singularity} containers/assembly-by-refseq.sif ")
semibin = stringr::str_glue("{singularity} run -B {basedir} containers/assembly-by-refseq.sif ")

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
variables.pigz = stringr::str_glue(" -p {threads} ")
variables.megahit = stringr::str_glue(" -t {threads}  --memory 0.8 ")
variables.flash = stringr::str_glue(" -t 4 -M 150 ")
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
  bbtools = "",
  pigz = ""
)

# DATABASE REFERENCES -----------------------------------------------------

# DECORATIONS
indent1 = "#---- "
indent2 = "#-------- "
indent3 = "#------------ "
indent4 = "#---------------- "
