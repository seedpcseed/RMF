#' RMetaflow
#' Module: fastp
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
  printGlue("#---- Scrub-a-dub-dub, some fastq in a tub", log = TRUE)
  
  # set up step variables ---------------------------------------------------
  
  printGlue("#-------- Getting directory information for {names(commands[STEP][1])}", log = TRUE)
  in.dir = dirNormalize(current.step[[1]]$from)
  out.dir = dirNormalize(names(current.step[1]))
  dirIf(out.dir)
  
  ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE) %>% sort()
  out.samples = sapply(list.files(out.dir, pattern = current.step[[1]]$out_ext), 
                       function(x) str_split(x, split_ext)[[1]][1]) %>% as.vector()
  
  if( length(out.samples) > 0 ){
    ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
    remaining.samples = samples[!samples%in%out.samples]
  }else{
    remaining.samples = samples
  }
  
  if( length(ls) > 0 ){
    printGlue("#-------- Preparing samples for fastp cleanup.", log = TRUE) 
    fseq = ls[ seq(from = 1, to = length(ls), by = 2)]
    rseq = ls[ seq(from = 2, to = length(ls), by = 2)]
    fout = paste0(out.dir, "/", remaining.samples, "_R1.fastq.gz")
    rout = paste0(out.dir, "/", remaining.samples, "_R2.fastq.gz")
    hout =  paste0(out.dir, "/", remaining.samples, ".html")
    jout =  paste0(out.dir, "/", remaining.samples, ".json")
    
    job.list = list(fseq = fseq, rseq = rseq,
                    fout = fout, rout = rout, 
                    hout = hout, jout = jout)
    
    calls <- purrr::pmap(job.list, function(fseq = fseq, rseq = rseq, fout =fout, rout = rout, hout = hout, jout = jout){
      str_glue("{fastp} fastp {variables.fastp} \\
                -i {fseq} \\
                -I {rseq} \\
                -o {fout} \\
                -O {rout} \\
                --detect_adapter_for_pe \\
                -j {jout} \\
                -h {hout} 2>&1 | tee  -a log/fastp.log")
                          
      })
    
    if( slurm ){
      printGlue("#-------- Oh I love arraying on your Slurm", log = TRUE)
      fastp.slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},  
      plan = "submit", 
      overwrite = TRUE, 
      njobs = length(calls),
      mc.cores = 16,
      sbatch_opt = list(partition = partition, account = account, time = "06:00:00"),
      preamble = c("module load singularity", "module load R/4.0.3"),
      job_name = paste0("Slurm-", projectname, "-fastp"), 
      tmp_path = projectdir)
      
      wait_slurm(fastp.slurm, freq = 60, timeout = -1)
      
      Slurm_clean(fastp.slurm)
    }else{
      printGlue("#--------- Sometime you just have to go one at a time (no Slurm)", log =TRUE)
      lapply(calls, function(x) {print(x); system(x)})  
    }
  }
  
  printGlue("#-------- Yes! I fastp'd your fastqs just for you!", log = TRUE)
  ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE) %>% sort()
  out.samples = sapply(list.files(out.dir, pattern = ".fastq.gz"), 
                       function(x) str_split(x, "_R")[[1]][1]) %>% as.vector()
  
  if( length(out.samples) > 0 ){
    ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
    if( length(ls) == 0 ){
      commands[[STEP]]$completed = TRUE
      updateYaml()
    }
  }
  
}


source.complete = TRUE