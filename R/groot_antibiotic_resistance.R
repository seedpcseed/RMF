#' RMetaflow
#' Module: humann_pathways
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
  printGlue("#---- I am groot", log = TRUE)

  # set up step variables ---------------------------------------------------
  
  in.dir = dirNormalize(current.step[[1]]$from)
  out.dir = dirNormalize(names(current.step[1]))
  dirIf(out.dir)
  
# find files and send command ---------------------------------
  
  ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE, recursive = TRUE) %>%
    str_subset(pattern = "intermediate", negate = TRUE) %>%
    sort()
  
  remaining.samples = samples

  printGlue("#-------- Prepare samples for groot analysis", log = TRUE)
  
  job.list = list(ls = ls, remaining.samples = remaining.samples)
  
  calls = purrr::pmap(job.list, function(ls, remaining.samples){
    str_glue("gunzip -c {ls} | \\
             {groot} groot align {variables.groot} > \\
             {out.dir}/{remaining.samples}.groot.bam && \\
             {groot} groot report -p {threads} --bamFile {out.dir}/{remaining.samples}.groot.bam > \\
             {out.dir}/remaining.samples.groot.report.txt -c 0.97")
  })
   
  if( slurm ){

  groot_slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},
                                               plan = "none",
                                               overwrite = TRUE,
                                               mc.cores = threads,
                                               njobs = length(calls),
                                               sbatch_opt = list(partition = partition, account = account, time = "08:00:00"),
                                               preamble = c("module load singularity", "module load R/4.0.3"),
                                               job_name = paste0("Slurm-", projectname, "-groot"),
                                               tmp_path = projectdir)

    sbatch(groot_slurm)

    wait_slurm(groot_slurm, freq = 60, timeout = -1)

    Slurm_clean(groot_slurm)

  }else{
    lapply(calls, function(x) {print(x[[1]]); system(x[[1]])})
  }
}

    ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE) %>% sort()
    out.samples = sapply(list.files(out.dir, pattern = current.step[[1]]$out_ext, recursive = TRUE),
                         function(x) str_split(x, ".contigs")[[1]][1]) %>% as.vector()

    if( length(out.samples) > 0 ){
      ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
      if( length(ls) ==0 ){
        commands[STEP][[1]]$completed = TRUE
        updateYaml()
      }

}


source.complete = TRUE