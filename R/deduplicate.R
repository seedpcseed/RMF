#' RMetaflow
#' Module: deduplicate
#' Initiated March 2021
#' Patrick C. Seed
#' Stanley Manne Children's Research Institute
#' Ann and Robert H. Lurie Children's Hospital
#' Northwestern University
#' MIT license
#' 

if ( !commands[STEP][[1]]$completed ){
  
  library(tidyverse)
  library(slurmR)
  
  printGlue("#---- Running {names(commands[STEP][1])}", log = TRUE)
  printGlue("#---- You are going to want to stick around for this one", log = TRUE)

  # set up step variables ---------------------------------------------------
  
  in.dir = dirNormalize(current.step[[1]]$from)
  out.dir = dirNormalize(names(current.step[1]))
  dirIf(out.dir)
  
# find files and send deduplicate command ---------------------------------

ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE) %>% sort()
out.samples = sapply(list.files(out.dir, pattern = current.step[[1]]$out_ext), 
                     function(x) str_split(x, "_R")[[1]][1]) %>% as.vector()

if( length(out.samples) > 0 ){
  ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
  remaining.samples = samples[!samples%in%out.samples]
}else{
  remaining.samples = samples
}

if( length(ls) > 0 ){
  
  fseq = ls[ seq(from = 1, to = length(ls), by = 2)]
  rseq = ls[ seq(from = 2, to = length(ls), by = 2)]
  temp.dir = paste0(out.dir, "/", remaining.samples)
  fseq.reloc = paste0(out.dir, "/", remaining.samples, "/", remaining.samples, "_R1.fastq")
  rseq.reloc = paste0(out.dir, "/", remaining.samples, "/", remaining.samples, "_R2.fastq")
  fout = paste0(out.dir, "/", remaining.samples, "_R1.fastq")
  rout = paste0(out.dir, "/", remaining.samples, "_R2.fastq")
  
  
  job.list = list(fseq = fseq, rseq = rseq,
                  temp.dir = temp.dir, 
                  fseq.reloc = fseq.reloc, rseq.reloc = rseq.reloc,
                  fout = fout, rout = rout)
  
  printGlue("#-------- Preparing samples for deduplication \n
            #-------- (Did I mention I was preparing samples for deduplication) \n
            #-------- Just want to remind you that I'm removing duplicates and \n
            #-------- I won't say it again", log = TRUE)
  calls <- purrr::pmap(job.list, function(fseq, rseq,temp.dir, fseq.reloc, rseq.reloc, fout, rout){
    str_glue("mkdir -p {temp.dir} &&  echo 'Opened directory {temp.dir}' && \\
    {pigz} pigz -dc {fseq} > {fseq.reloc} && echo 'Decompressed {fseq.reloc}' && \\
    {pigz} pigz -dc {rseq} > {rseq.reloc} && echo 'Decompressed {rseq.reloc}' && \\
    echo 'Starting deduplication' && {idseq.dedup} idseq-dedup {variables.idseq.dedup} -i {fseq.reloc} -i {rseq.reloc} \\
    -o {fout} -o {rout} 2>&1 | tee -a log/deduplicate.log && echo 'Completed deduplication' && \\
    echo 'Compressing {fout} and {rout}' && {pigz} pigz {fout} {rout} && \\
    rm -r {temp.dir} && echo 'Completed deduplication' ")})

  printGlue("Sending calls in {names(commands[STEP][1])}", log = TRUE)
  
  if( slurm ){
    printGlue("Where did the term Slurm come from? Anyway, I'm going with it", log = TRUE)
    dedup.slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},  
    plan = "none", 
    overwrite = TRUE, 
    mc.cores = 1,
    njobs = length(calls),
    sbatch_opt = list(partition = partition, account = account, time = "06:00:00", mem = "8G"),
    preamble = c("module load singularity", "module load R/4.0.3"),
    job_name = paste0("Slurm-", projectname, "-dedup"), 
    tmp_path = projectdir)
    
    sbatch(dedup.slurm)
    
    wait_slurm(dedup.slurm, freq = 60, timeout = -1)
    
    Slurm_clean(dedup.slurm)
  
    }else{
      printGlue("#-------- Deduplicating one file at a time. Sorry about that.....")
      lapply(calls, function(x) {print(x); system(x)})  
    }
}

  printGlue("#-------- I've deduplicated them all for you and now I don't have to \n
            #-------- repeat myself", log = TRUE)
  ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE) %>% sort()
  out.samples = sapply(list.files(out.dir, pattern = current.step[[1]]$out_ext), 
                       function(x) str_split(x, "_R")[[1]][1]) %>% as.vector()
  
  if( length(out.samples) > 0 ){
    ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
    if( length(ls) ==0 ){
      commands[[STEP]]$completed = TRUE
      updateYaml()
    }

  }
}



source.complete = TRUE