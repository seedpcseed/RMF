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
  printGlue("#---- Finding the microbe in your humann")

  # set up step variables ---------------------------------------------------

  in.dir = dirNormalize(current.step[[1]]$from)
  out.dir = dirNormalize(names(current.step[1]))
  dirIf(out.dir)

  # find files and send command ---------------------------------

    ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE, recursive = TRUE) %>%
      str_subset(pattern = "intermediate", negate = TRUE) %>%
      sort()

    remaining.samples = samples

    printGlue("#------- Prepare samples for humann pathway analysis", log = TRUE)
    printGlue("#------- Some humans are slow. This humann will be quick", log = TRUE)

    jobs.list = list(fseq = ls[seq(1, length(ls), by = 2)], rseq = ls[seq(2, length(ls), by = 2)], remaining.samples = samples)
    calls = purrr::pmap(jobs.list, function(fseq, rseq, remaining.samples){
      str_glue("cat {fseq} {rseq} > {out.dir}/{remaining.samples}.concat.fastq.gz && \\
               {classify} humann {variables.humann} --input {out.dir}/{remaining.samples}.concat.fastq.gz --output {out.dir}/{remaining.samples}")
    })

    if( slurm ){
    printGlue("#-------- Did you ever know Slurm's your hero? ", log = TRUE)
    humann_slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},
                                                 plan = "none",
                                                 overwrite = TRUE,
                                                 mc.cores = threads,
                                                 njobs = length(calls),
                                                 sbatch_opt = list(partition = partition, account = account, time = "24:00:00",
                                                                   'mem-per-cpu' = "4G", 'cpus-per-task' = threads),
                                                 preamble = c("module load singularity", "module load R/4.0.3"),
                                                 job_name = paste0("Slurm-", projectname, "-humann"),
                                                 tmp_path = projectdir)

      sbatch(humann_slurm)

      wait_slurm(humann_slurm, freq = 60, timeout = -1)

      Slurm_clean(humann_slurm)

    }else{
      printGlue("#-------- We'll move through your samples in serial but \n
                your samples won't get soggy", log = TRUE)
      lapply(calls, function(x) {print(x[[1]]); system(x[[1]])})
    }


  # cleanup files -----------------------------------------------------------
  ls = list.files(out.dir, pattern = "concat", full.names = TRUE)
  lapply(ls, file.remove)

  ls = list.files(out.dir, all.files = TRUE, recursive = TRUE, full.names = TRUE)
  lapply(ls, function(x) file.copy(from = x, to = out.dir, copy.mode = FALSE))


  ld = list.dirs(out.dir, full.names = TRUE)[-1]
  lapply(ld, function(x) system(str_glue("rm -r {x}")))


  # merge tables for many samples -------------------------------------------

  calls = list(str_glue("{classify} humann_join_tables -i {out.dir} -o {out.dir}/merged.genefamiles.tsv \\
                        --file_name 'gene' && \\
                        {classify} humann_renorm_table --input {out.dir}/merged.genefamiles.tsv \\
                        --units 'cpm' --output {out.dir}/merged.genefamiles.cpm_norm.tsv"),
               str_glue("{classify} humann_join_tables -i {out.dir} -o {out.dir}/merged.pathcoverage.tsv \\
                        --file_name 'coverage'"),
               str_glue("{classify} humann_join_tables -i {out.dir} -o {out.dir}/merged.pathabundance.tsv \\
                        --file_name 'abundance' && \\
                        {classify} humann_renorm_table --input {out.dir}/merged.pathabundance.tsv \\
                        --units 'cpm' --output {out.dir}/merged.pathabundance.cpm_norm.tsv"))

  if( slurm ){
    printGlue("#-------- Did you ever know Slurm's your hero? ", log = TRUE)
    humann_slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},
                                 plan = "none",
                                 overwrite = TRUE,
                                 mc.cores = 1,
                                 njobs = length(calls),
                                 sbatch_opt = list(partition = partition, account = account, time = "24:00:00"),
                                 preamble = c("module load singularity", "module load R/4.0.3"),
                                 job_name = paste0("Slurm-", projectname, "-humann"),
                                 tmp_path = projectdir)

    sbatch(humann_slurm)

    wait_slurm(humann_slurm, freq = 60, timeout = -1)

    Slurm_clean(humann_slurm)

  }else{
    printGlue("#-------- We'll move through your samples in serial but \n
                your samples won't get soggy", log = TRUE)
    lapply(calls, function(x) {print(x[[1]]); system(x[[1]])})
  }

  # clean up intermediate single sample files -------------------------------
  ls = list.files(out.dir, pattern = "concat", full.names = TRUE)
  dir.create(file.path(out.dir, "single_sample_files"), showWarnings = FALSE)
  lapply(ls, function(x) system(str_glue("mv {x} {out.dir}/single_sample_files")))

  ls = list.files(file.path(out.dir, "single_sample_files"), full.names = TRUE)
  calls = str_glue("{classify} pigz {ls}")

  if( slurm ){
    printGlue("#-------- you should send Slurm a thank you card \n
              the way it's workin' for you.....", log = TRUE)
    humann_slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},
                                 plan = "none",
                                 overwrite = TRUE,
                                 mc.cores = 1,
                                 njobs = length(calls),
                                 sbatch_opt = list(partition = partition, account = account, time = "24:00:00"),
                                 preamble = c("module load singularity", "module load R/4.0.3"),
                                 job_name = paste0("Slurm-", projectname, "-humann"),
                                 tmp_path = projectdir)

    sbatch(humann_slurm)

    wait_slurm(humann_slurm, freq = 60, timeout = -1)

    Slurm_clean(humann_slurm)

  }else{
    printGlue("#-------- Zippidee do dah, zippidee eh \n
              zipping one file at a time might take much of the \n
              day ......", log = TRUE)
    lapply(calls, function(x) {print(x[[1]]); system(x[[1]])})
  }

  # check if done with this module ------------------------------------------

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

}

source.complete = TRUE
