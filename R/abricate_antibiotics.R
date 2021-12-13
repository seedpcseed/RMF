#' RMetaflow
#' Module: abricate antibiotic resistome
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
  library(R.utils)
  
  printGlue("#---- Running {names(commands[STEP][1])}", log = TRUE)
  printGlue("#---- You cannot resist!")

  # set up step variables ---------------------------------------------------
  
  in.dir = dirNormalize(current.step[[1]]$from)
  out.dir = dirNormalize(names(current.step[1]))
  dirIf(out.dir)
  
  # find files and send command ---------------------------------
    
    ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE, recursive = TRUE) %>%
      str_subset(pattern = "intermediate", negate = TRUE) %>%
      sort()
    
    remaining.samples = samples
  
    printGlue("#------- Let's find some resistance, shall we?", log = TRUE)
    
    calls = lapply(ls, function(x){
      sample.string = str_split(basename(x), ".contigs")[[1]][1]
      str_glue("{abricate} abricate {variables.abricate} {x} > {out.dir}/{sample.string}.resistance.abricate.tsv && \\
               {abricate} abricate --threads {threads} --db vfdb --nopath {x} > {out.dir}/{sample.string}.vf.abricate.tsv ") %>% 
        str_squish()
    })
    
    
    if( slurm ){
    printGlue("#-------- Did you ever know Slurm's your hero? ", log = TRUE)  
    abricate_slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},  
                                                 plan = "none", 
                                                 overwrite = TRUE, 
                                                 mc.cores = threads,
                                                 njobs = length(calls),
                                                 sbatch_opt = list(partition = partition, account = account, time = "24:00:00", mem = "64G"),
                                                 preamble = c("module load singularity", "module load R/4.0.3"),
                                                 job_name = paste0("Slurm-", projectname, "-abricate"), 
                                                 tmp_path = projectdir)
      
      sbatch(abricate_slurm)
      
      wait_slurm(abricate_slurm, freq = 60, timeout = -1)
      
      Slurm_clean(abricate_slurm)
      
    }else{
      printGlue("#-------- We'll move through your samples in serial but \n
                your samples won't get soggy", log = TRUE)
      lapply(calls, function(x) {print(x[[1]]); system(x[[1]])})  
    }


  # make a summary table for all samples ------------------------------------
  ls = list.files(out.dir, pattern = "resistance.abricate.tsv", full.names = TRUE)
  calls = str_glue("{abricate} abricate --summary {paste(ls, collapse = ' ')} > {out.dir}/all.samples.resistance.abricate.tab ")
  system(calls)
  
  ls = list.files(out.dir, pattern = "vf.abricate.tsv", full.names = TRUE)
  calls = str_glue("{abricate} abricate --summary {paste(ls, collapse = ' ')} > {out.dir}/all.samples.vf.abricate.tab ")
  system(calls)
      
  # cleanup files -----------------------------------------------------------
  
  # check if done with this module ------------------------------------------
  
  if( length(out.samples) > 0 ){
    ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
    if( length(ls) ==0 ){
      commands[STEP][[1]]$completed = TRUE
      updateYaml()
    }
  
  }

}

source.complete = TRUE