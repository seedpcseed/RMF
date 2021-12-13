#' RMetaflow
#' Module: classify (kaiju)
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
  printGlue("#---- Nice taxonomy, My Friend!")

  # set up step variables ---------------------------------------------------

  in.dir = dirNormalize(current.step[[1]]$from)
  out.dir = dirNormalize(names(current.step[1]))
  dirIf(out.dir)

# find files ---------------------------------

ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE) %>% sort()
out.samples = sapply(list.files(out.dir, pattern = current.step[[1]]$out_ext),
                     function(x) str_split(x, "_kaiju.tsv")[[1]][1]) %>% as.vector()

if( length(out.samples) > 0 ){
  ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
  remaining.samples = samples[!samples%in%out.samples]
}else{
  remaining.samples = samples
}

if( length(ls) >0 ){

  fseq = ls[ seq(from = 1, to = length(ls), by = 2)]
  rseq = ls[ seq(from = 2, to = length(ls), by = 2)]
  kout = paste0(out.dir, "/", remaining.samples, "_kaiju.out")
  taxout = paste0(out.dir, "/", remaining.samples, "_kaiju.taxa")
  taxtable = paste0(out.dir, "/", remaining.samples, "_kaiju.tsv")


  job.list = list(fseq = fseq, rseq = rseq,
                  kout = kout, taxout =  taxout, taxtable = taxtable)

  printGlue("#-------- Preparing samples for kaiju classification", log = TRUE)

  calls <- purrr::pmap(job.list, function(fseq, rseq, kout, taxout, taxtable){
    str_glue("{classify} kaiju-multi {variables.kaiju.map} -i {fseq} -j {rseq} -o {kout} 2>&1 | \\
    tee -a log/classify.log && \\
    {classify} kaiju-addTaxonNames {variables.kaiju.tax} -i {kout} \\
    -o {taxout} 2>&1 | tee  -a log/classify.log && \\
    {classify} kaiju2table {variables.kaiju.table} \\
    -o {taxtable} {taxout} 2>&1 | tee -a log/classify.log")})

  if( slurm ){
    printGlue("#-------- Sending those babies to Slurm and churn", log = TRUE)
    classify.slurm <- Slurm_lapply(calls, function(x) {print(x[[1]]); system(x[[1]])},
    plan = "none",
    overwrite = TRUE,
    mc.cores = threads,
    njobs = length(calls),
    sbatch_opt = list(partition = partition, account = account, time = "04:00:00"),
    preamble = c("module load singularity", "module load R/4.0.3"),
    job_name = paste0("Slurm-", projectname, "-classify"),
    tmp_path = projectdir)

    sbatch(classify.slurm)

    wait_slurm(classify.slurm, freq = 60, timeout = -1)

    Slurm_clean(classify.slurm)

    }else{
      printGlue("#-------- No Slurm here but still will churn (just one after another)", log = TRUE)
      lapply(calls, function(x) {print(x); system(x)})
    }

    gather_kaiju_tables(in.dir = out.dir)

    call = str_glue("{pigz} pigz {variables.pigz} {out.dir}/*_kaiju.out {out.dir}/*_kaiju.taxa")
    print(call)
    system(call)
}

  printGlue("#--------- Whew! Naming all of those precious babies was a lot of work", log = TRUE)
  printGlue("#--------- but on to our next job!", log = TRUE)
  printGlue("#------------- (pst... are you having fun yet?")
  ls = list.files(in.dir, pattern = current.step[[1]]$in_ext, full.names = TRUE) %>% sort()
  out.samples = sapply(list.files(out.dir, pattern = current.step[[1]]$out_ext),
                       function(x) str_split(x, "_kaiju.tsv")[[1]][1]) %>% as.vector()

  if( length(out.samples) > 0 ){
    ls = ls[!grepl(paste0(out.samples, collapse = "|"), ls)]
    if( length(ls) ==0 ){
      commands[[STEP]]$completed = TRUE
      updateYaml()
    }

  }
}



source.complete = TRUE
