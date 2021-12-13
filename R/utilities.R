#' RMetaflow
#' Module: utilities.R
#' Initiated March 2021
#' Patrick C. Seed
#' Stanley Manne Children's Research Institute
#' Ann and Robert H. Lurie Children's Hospital
#' Northwestern University
#' MIT license
#' 

library(tidyverse)
library(rprojroot)
library(yaml)

# general utilities -------------------------------------------------------

getProjDir = function(sentinel = "startup.sh"){
  return(find_root(has_file(sentinel)))
}

getCurrentStep = function(){
  for(i in 1:length(commands)){
    if( !commands[[i]]$completed) {
      STEP <<- i 
      return("Loaded steps")
    } 
  }
}

printGlue <- function(in.string, log = FALSE){
  print(str_glue(in.string)); cat("\n")
  if(log){
    st = as.character(Sys.time())
    converted.string = toString(str_glue(str_glue("{st}] {in.string}")))
    gsub("\r?\n|\r", " ", str_glue("echo '{converted.string}' 2>&1 | tee -a log/log.txt")) %>% system() 
  }
}

updateYaml <- function(){
  write_yaml(commands, file = "yaml/metaflow-steps.yml")
}

dirNormalize <- function(in.string){
  if( length(str_split(in.string, " ")[[1]])> 1){in.string = str_split(in.string, " ")[[1]][1]}
  return(paste0(str_split(in.string, "\\.")[[1]], collapse = "/"))
}
  
stepIf <- function(in.source){
  if( length(str_split(in.source, ":::")[[1]]) != 2 ){
    str_glue("{Sys.time()}] Command does not have the necessary 2 parts:
          TYPE ::: COMMAND")
    q('no')
  }
  current.step$type <<- str_split(in.source, ":::")[[1]][1]
  current.step$call <<- str_split(in.source, ":::")[[1]][2]
}

dirIf <- function(in.dir){
  if(!dir.exists(in.dir)){
    dir.create(in.dir, recursive = TRUE)
    print(str_glue("Directory created: {in.dir}"))
  }
}

fileIf <- function(in.command, in.exists = FALSE){
  if( length(str_split(in.command, ":::")[[1]]) != 2){
    str_glue("{Sys.time()}] Command does not have the necessary 2 parts:
          TYPE ::: COMMAND")
    q('no')
  }
  command.type = str_split(in.command, ":::")[[1]][1]
  command.to.issue = str_split(in.command, ":::")[[1]][2]
  if ( command.type == "system"){
    system(toString(command.to.issue))
  }
  if ( command.type == "source"){
    source(command.to.issue)
  }
  printGlue(command.to.issue, log = TRUE)
}

markCompleted <- function(completed.step){
  commands[STEP][[1]]$completed = TRUE
}

getCommandSteps <- function(yaml.in = "yaml/metaflow-steps.yml"){
  commands <<- yaml.load_file(yaml.in)
  names(commands)
}

findPatterns <- function(in.dir){
  # script will return 3 global variables: best.split, pairs, samples
  # based on best guess of sample name nomenclature
  # first check for pairs
  ls = list.files(in.dir)
  if (length(ls) %% 2 != 0 ){
    printGlue("Input directory {in.dir} does not have an even number of files.
              Have to assume there is something wrong or the files are not pair-end.
              Please check to continue.")
    q('no')
  }
  
  lsR = list.files(in.dir, pattern = "_R")
  lsL = list.files(in.dir, pattern = "_L")
  ls1 = list.files(in.dir, pattern = "_1")
  
  max.length = max(length(lsR), length(lsL), length(ls1))
  max.position = which(c(length(lsR), length(lsL), length(ls1))==max.length)
  
  best.split <<- case_when(max.position == 1 ~ "_R",
            max.position == 2 ~ "_L",
            max.position == 3 ~ "_1")
  
  if ( best.split == "_1" ){
    pair_ext <<- c("_1", "_2")
  }else{
    pair_ext <<- c("_R1", "_R2")
  }
  
  ls = list.files(in.dir, pattern = best.split)
  samples <<- sapply(ls, function(x) str_split(x, best.split)[[1]][1]) %>% as.vector() %>% unique()
  
}

gather_kaiju_tables<- function(in.dir, out.name = "combined.kaiju", in.pattern = "tsv"){
  in.files = list.files(in.dir, pattern = in.pattern, full.names = TRUE) %>%
    str_subset(pattern = "combined", negate = TRUE)
  df1 = read_tsv(in.files[1], col_types = cols())  %>% filter(!is.na(taxon_id)) %>% as.data.frame()
  df = data.frame(matrix(nrow = 0 , ncol=length(df1$taxon_id)), check.names = FALSE)
  colnames(df) = df1$taxon_id
  
  taxa = data.frame(matrix(nrow = 0, ncol = 8))
  colnames(taxa) = c("taxon_id", "Kingdom", "Phylum", "Class", 
                     "Order", "Family", "Genus", "Species")
  
  for(x in 1:length(in.files)){
    svMisc::progress((x-1), max.value = (length(in.files)-1) , progress.bar = TRUE)
    df.temp <- data.table::fread(in.files[x]) %>% as.data.frame() %>% filter(!is.na(taxon_id))
    rownames(df.temp) = df.temp$taxon_id
    sample = basename(df.temp$file[1])
    sample.hash  = digest::digest(sample, "crc32")
    df.sub = df.temp %>% select(reads)
    row.temp = data.frame(t(df.sub), row.names = str_glue("{sample}"), check.names = FALSE) 
    row.temp = row.temp[complete.cases(row.temp)]
    df = plyr::rbind.fill(df, row.temp)
    
    taxa.temp = df.temp %>% select(taxon_id, taxon_name) %>% separate(taxon_name, c("Superkingdom", "Phylum", "Class", 
                                                                                    "Order", "Family", "Genus", "Species"),
                                                                      sep = ";")
    taxa.temp <- taxa.temp[!duplicated(taxa.temp$taxon_id, taxa$taxon_id),]
    
    taxa = rbind(taxa, taxa.temp)
  }
  
  rownames(df) = sapply(basename(in.files), function(x) str_split(x, ".kaiju")[[1]][1])
  df = data.frame(Sample = rownames(df), df, check.names = FALSE)
  df[is.na(df)] <- 0
  
  taxa <- taxa %>% distinct()
  write.table(df, file=str_glue("{in.dir}/{out.name}.counts.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(taxa, file=str_glue("{in.dir}/{out.name}.taxa.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
}

fix.orna <- function(infile, suffix = ".normalized.fa"){
  count.lines = countLines(infile)[1]
  outfile = paste0(str_split(infile, suffix)[[1]][1], ".normalized.corrected.fa")
  for(i in seq(from = 0, to = (count.lines-1), by = 1)){
    inlines = read_lines(opt$infile, skip = i, n_max = 1) %>% str_replace_all("@", ">") %>% str_replace_all("\\+", "")
    if( nchar(inlines) != 0){
      write_lines(inlines, file = outfile, append = TRUE)
    }
  }
}

reset.startup <- function(in.step){
  source("R/utilities.R")
  STEP <<- in.step
  current.step <<- commands[in.step]
  dirIf(dirNormalize(names(current.step)))
  stepIf(current.step[[1]]$source)
  
  in.dir <<- dirNormalize(current.step[[1]]$from)
  out.dir <<- dirNormalize(names(current.step[1]))
  dirIf(out.dir)
  
  printGlue("# --------- COMMAND ORDER ")
  print(names(commands))
  printGlue("# ----------------------------------")
  printGlue("# --------- in.dir = {in.dir}")
  printGlue("# --------- out.dir = {out.dir}")
  printGlue("# --------- current.step ")
  print(unlist(current.step))
}
  
na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

add_slurm_throttle <- function(in.file, throttle = 1){
  file.read = readr::read_lines(in.file)
  array.line = file.read[grep("array", file.read)]
  new.array.line = paste0(array.line, str_glue("%{throttle}"))
  file.read[grep("array", file.read)] <- new.array.line
  write_lines(file.read, in.file)
}


findPackages <- function(file){
  txt <- readLines(file)
  inx <- grep('::', txt)
  txt <- txt[inx]
  m <- regexpr('[[:alnum:]]+::', txt)
  pkg <- regmatches(txt, m)
  unique(sub('::', '', pkg))
}

find_needed_packages <- function(){
  pkgs.to.install <- attachment::att_from_rscripts("R")
  install.packages(pkgs.to.install, dependencies = TRUE)
}

split_contigs <- function(infile, outdir, clines = 10000, t = 8){
  
  line.marker = 0
  total.lines = system(str_glue("cat {infile} | wc -l"), intern = TRUE) %>% as.numeric()
  system(str_glue("touch {outdir}/contigs.separated.fasta"))
  
  while((line.marker + clines) < total.lines){
    
    if((line.marker + clines) < total.lines){
      add.lines = clines
    }else{
        add.lines = total.lines - line.marker
    }
    
    print(str_glue("Current Lines {line.marker} - {line.marker + add.lines}"))
    infile.read <- read_lines(infile, n_max = (add.lines + line.marker), skip = line.marker, num_threads = t, lazy = FALSE)
    contig.start <- grep(">", infile.read)
    if(length(contig.start>0)){infile.read <- infile.read[-contig.start]}
    infile.read <- c(str_glue("> contig_{line.marker}"), infile.read)
    write_lines(infile.read, str_glue("{outdir}/contigs.separated.fasta"), append = TRUE)
    line.marker = line.marker + clines + 1
  }
  
}

split.fasta <- function(infile, outdir, unit.size = 10000, t = 8){
  total.lines = system(str_glue("cat {infile} | wc -l"), intern = TRUE) %>% as.numeric()
  n.parts = floor(total.lines/unit.size)
  line.insert.points <- seq(from = 1, to = total.lines, by = unit.size)
  system(str_glue("cp {infile} {outdir}/split.{unit.size}.fna"))
  sapply(line.insert.points, function(x){
    call = str_glue("sed -i '{x}i >contig_{x}' {outdir}/split.{unit.size}.fna")
    system(call)
    Sys.sleep(3)
    
  })
  
  #system(str_glue("split -l {n.parts} {infile} {outdir}/parts"))
  #ls = list.files(str_glue("{outdir}"), pattern = "parts", full.names = TRUE)
  #system(str_glue("touch {outdir}/contigs.separated.fasta"))
  
  for(i in seq_along(ls)){
    infile.read <- read_lines(ls[i], num_threads = t, lazy = FALSE, progress = FALSE)
    contig.start <- grep(">", infile.read)
    if(length(contig.start>0)){infile.read <- infile.read[-contig.start]}
    infile.read <- c(str_glue("> contig_{i}"), infile.read)
    write_lines(infile.read, str_glue("{outdir}/contigs.separated.fasta"), append = TRUE)
    print(str_glue("Completed: {ls[i]}"))
    file.remove(ls[i])
  }

    system(str_glue("rm {outdir}/parts*"))
}

strip.file.PreSuffixes <-function(in.file){
  extensions <- ".fasta|.fa|.fna|.txt|.gz|.fq|.fastq| 
                  .faa|.tar|.tsv|.csv|.out" %>% stringr::str_squish()
  stripped1 = tools::file_path_sans_ext(basename(in.file))
  stripped.check = stringr::str_detect(stripped1, extensions)
  if(stripped.check){stripped1 = tools::file_path_sans_ext(basename(stripped1))}
  return(stripped1)
}

make_custom_referenceseeker_db <- function(db_name, db_dir, genome_dir, taxa = "unspecified", tax_id, 
                                           genome_ext = "fna.gz", genome_status = "complete"){
  #initiate the database
  system(str_glue("{assembly.by.refseq} referenceseeker_db init --output {db_dir} \\
  -d {db_name}") %>% str_squish())
  
  #relocate genomes if needed
  if( genome_dir != str_glue("{db_dir}/{db_name}")){
    all.genomes = list.files(str_glue("{genome_dir}"), pattern = genome_ext, full.names = TRUE, recursive = TRUE)
    calls = str_glue("mv {all.genomes} {db_dir}/{db_name}")
    future.apply::future_lapply(calls, system)
  }
  
  all.genomes = list.files(str_glue("{db_dir}/{db_name}"), pattern = genome_ext)
  
  system(str_glue("[ ! -e {db_dir}/assembly-by-refseq.sif ] && cp containers/assembly-by-refseq.sif {db_dir}/assembly-by-refseq.sif"))
  
  lapply(all.genomes, function(x){
    calls = str_glue("cd {db_dir} && gzip -d {db_name}/{x}; \\
    singularity run -B {basedir} assembly-by-refseq.sif referenceseeker_db import --db {db_name} --genome {db_name}/{x %>% str_replace(genome_ext, 'fna')} \\
             --taxonomy {tax_id} --organism {taxa} --status {genome_status}")
    system(calls)
  })
  
  system(str_glue("rm {db_dir}/assembly-by-refseq.sif"))
  
}

get_ncbi_genomes <- function(genus_name, genome_type = "genbank", genome_format = "fasta",
                             genome_group = "bacteria", output.dir){
  calls = str_glue("{assembly.by.refseq} ncbi-genome-download -s {genome_type} -F {genome_format} \\
                   --genera {genus_name} -o {output.dir} {genome_group}")
  system(calls)
}


